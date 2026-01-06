# Tests for LaTeX output functionality

@testset "LaTeX Type Construction" begin
    @variables a b c
    
    # Test LaTeX wrapper construction with different types
    @test LaTeX(a) isa LaTeX
    @test LaTeX([a, b, c]) isa LaTeX
    @test LaTeX([a b; b c]) isa LaTeX
    @test LaTeX(a + b) isa LaTeX
    @test LaTeX(sqrt(a)) isa LaTeX
    
    # Test that expr field is stored correctly
    l = LaTeX(a)
    @test l.expr === a
    
    lv = LaTeX([a, b])
    @test isequal(lv.expr, [a, b])
    
    M = [a b; b c]
    lm = LaTeX(M)
    @test isequal(lm.expr, M)
end

@testset "LaTeX Plain Text Display" begin
    @variables x y
    
    # Test plain text fallback display
    l = LaTeX(x + y)
    io = IOBuffer()
    show(io, MIME("text/plain"), l)
    output = String(take!(io))
    @test !isempty(output)
    
    # Test regular show
    io = IOBuffer()
    show(io, l)
    output = String(take!(io))
    @test !isempty(output)
end

@testset "LaTeX HTML Display Without Latexify" begin
    @variables a b c
    
    # Reset the Latexify state to test fallback behavior
    # Note: We can't easily unload Latexify if it's loaded, so we test what we can
    
    # Test HTML display with vector
    lv = LaTeX([a, b, c])
    io = IOBuffer()
    show(io, MIME("text/html"), lv)
    html_output = String(take!(io))
    # Should contain some output (either LaTeX or fallback)
    @test !isempty(html_output)
    
    # Test HTML display with matrix
    M = [a b; b c]
    lm = LaTeX(M)
    io = IOBuffer()
    show(io, MIME("text/html"), lm)
    html_output = String(take!(io))
    @test !isempty(html_output)
    
    # Test HTML display with single expression
    l = LaTeX(a^2 + sqrt(b))
    io = IOBuffer()
    show(io, MIME("text/html"), l)
    html_output = String(take!(io))
    @test !isempty(html_output)
end

@testset "LaTeX MIME Display" begin
    @variables x y z
    
    # Test text/latex MIME type
    l = LaTeX(x + y)
    io = IOBuffer()
    show(io, MIME("text/latex"), l)
    latex_output = String(take!(io))
    @test !isempty(latex_output)
    
    # Test with vector
    lv = LaTeX([x, y, z])
    io = IOBuffer()
    show(io, MIME("text/latex"), lv)
    latex_output = String(take!(io))
    @test !isempty(latex_output)
    
    # Test with matrix
    M = [x y; y z]
    lm = LaTeX(M)
    io = IOBuffer()
    show(io, MIME("text/latex"), lm)
    latex_output = String(take!(io))
    @test !isempty(latex_output)
end

@testset "LaTeX with Eigenvalues" begin
    @variables a b
    
    # Create a symmetric matrix and get eigenvalues
    M = [a 0; 0 b]
    vals = eigvals(M)
    
    # Wrap in LaTeX
    l = LaTeX(vals)
    @test length(l.expr) == length(vals)
    
    # Test display
    io = IOBuffer()
    show(io, MIME("text/html"), l)
    html_output = String(take!(io))
    @test !isempty(html_output)
end

@testset "LaTeX Internal Functions" begin
    # Test _ensure_latexify function
    # This tests the lazy loading mechanism
    result = SymbolicDiagonalization._ensure_latexify()
    @test result isa Bool
    
    # Test _to_latex function
    @variables x
    latex_str = SymbolicDiagonalization._to_latex(x)
    @test latex_str isa String
    @test !isempty(latex_str)
    
    # Test with different expression types
    latex_str = SymbolicDiagonalization._to_latex(x^2 + 1)
    @test latex_str isa String
    
    # Test with numeric values
    latex_str = SymbolicDiagonalization._to_latex(42)
    @test latex_str == "42"
end
