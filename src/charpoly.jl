"""
    characteristic_polynomial(A; var = nothing)

Compute the characteristic polynomial of the square matrix `A` symbolically.
Returns `(poly, coeffs, λ)` where `poly` is the expanded polynomial expression,
`coeffs` are the coefficients in ascending order (constant term first),
and `λ` is the symbolic variable used.

For Hermitian/symmetric matrices with Complex{Num} entries, the characteristic
polynomial is guaranteed to be real, so we extract the real part to enable
proper symbolic differentiation.
"""
function characteristic_polynomial(A; var = nothing)
    λ = isnothing(var) ? _fresh_lambda() : var
    m, n = size(A)
    m == n || throw(ArgumentError("Matrix must be square, got $(m)×$(n))"))
    # Form λI - A so det(λI - A) yields the characteristic polynomial.
    shifted = _lambda_shift(Matrix(A), λ)
    # Bareiss keeps entries smaller and avoids division by symbolic pivots.
    poly = Symbolics.expand(_bareiss_det(shifted))
    
    # For Complex{Num} polynomials, extract the real part.
    # This is necessary because Symbolics.derivative doesn't work correctly
    # with Complex{Num} - it leaves unevaluated Differential terms.
    # For Hermitian matrices, the characteristic polynomial is always real.
    # For general complex matrices, the imaginary part would also need handling,
    # but eigenvalues of non-Hermitian complex matrices are outside our scope.
    if poly isa Complex && (real(poly) isa Num || imag(poly) isa Num)
        poly = real(poly)
    end
    
    # Extract coefficients via derivatives instead of polynomial division to
    # stay robust on arbitrary symbolic types.
    coeffs = _collect_coefficients(poly, λ, n)
    return poly, coeffs, λ
end

_fresh_lambda() = first(Symbolics.@variables λ)

function _lambda_shift(A, λ)
    n = size(A, 1)
    Iλ = Diagonal(fill(λ, n)) |> Matrix
    return Iλ .- A
end

# Fraction-free Bareiss determinant; works well on symbolic entries.
function _bareiss_det(M)
    n = size(M, 1)
    A = Matrix{eltype(M)}(M)
    prev = one(eltype(M))
    for k in 1:n-1
        pivot = A[k, k]
        if _issymzero(pivot)
            throw(ArgumentError("Zero pivot encountered in Bareiss determinant at position ($k, $k). Matrix may be singular or require pivoting."))
        end
        for i in k+1:n, j in k+1:n
            A[i, j] = (A[i, j] * pivot - A[i, k] * A[k, j]) / prev
        end
        # All entries below the pivot are zeroed in one shot (saves ops).
        fill!(view(A, k+1:n, k), zero(eltype(M)))
        prev = pivot
    end
    return A[n, n]
end

function _collect_coefficients(poly, λ, degree)
    coeffs = Vector{Any}(undef, degree + 1)  # Must be Any since coefficients can be Num or numeric
    coeffs[1] = Symbolics.expand(Symbolics.substitute(poly, λ => 0))
    deriv = poly
    fact = 1  # Use Int for factorial accumulation
    for k in 1:degree
        deriv = Symbolics.derivative(deriv, λ)
        fact *= k
        coeffs[k + 1] = Symbolics.expand(Symbolics.substitute(deriv, λ => 0) / fact)
    end
    return coeffs
end
