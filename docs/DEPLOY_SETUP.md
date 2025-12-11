# GitHub Pages Documentation Deployment Setup

Quick setup guide for deploying SymbolicDiagonalization.jl documentation to GitHub Pages.

## Prerequisites

- GitHub repository created
- Local repository with documentation files
- Julia installed locally

## Setup Steps

### Step 1: Generate SSH Keys

```bash
# Install DocumenterTools
julia -e 'using Pkg; Pkg.add("DocumenterTools")'

# Generate keys (run from repository root)
julia -e 'using DocumenterTools; DocumenterTools.genkeys()'
```

This will:
- Create `docs/.documenter` (private key)
- Display the public key in your terminal

**IMPORTANT**: Copy the public key - you'll need it in the next step!

### Step 2: Add Deploy Key to GitHub

1. Go to your GitHub repository
2. **Settings** → **Deploy keys** → **Add deploy key**
3. Title: `documenter-key`
4. Key: Paste the public key from Step 1
5. ✅ **Check "Allow write access"** (critical!)
6. Click **Add key**

### Step 3: Add Private Key to GitHub Secrets

1. Open `docs/.documenter` and copy its contents
2. Go to your repository: **Settings** → **Secrets and variables** → **Actions**
3. Click **New repository secret**
4. Name: `DOCUMENTER_KEY`
5. Value: Paste the entire contents of `docs/.documenter`
6. Click **Add secret**

### Step 4: Update Repository URLs

Edit `docs/make.jl` and replace `YOUR_USERNAME` with your GitHub username in TWO places:

```julia
# Line ~20
repo = "https://github.com/YOUR_USERNAME/SymbolicDiagonalization.jl/blob/{commit}{path}#{line}",
canonical = "https://YOUR_USERNAME.github.io/SymbolicDiagonalization.jl",

# Line ~30
deploydocs(;
    repo = "github.com/YOUR_USERNAME/SymbolicDiagonalization.jl",
```

### Step 5: Commit and Push

```bash
git add .github/workflows/documentation.yml
git add docs/make.jl
git add docs/.gitignore
git commit -m "Add GitHub Pages documentation deployment"
git push origin main
```

### Step 6: Enable GitHub Pages

1. Go to **Settings** → **Pages**
2. Source: **Deploy from a branch**
3. Branch: Select `gh-pages` (appears after first successful deployment)
4. Folder: `/ (root)`
5. Click **Save**

**Note**: The `gh-pages` branch won't exist until the first successful deployment runs.

### Step 7: Wait for First Deployment

1. Go to **Actions** tab
2. Watch the "Documentation" workflow run
3. Once complete (green checkmark), the `gh-pages` branch will exist
4. Go back to **Settings** → **Pages** and select `gh-pages` branch
5. Wait 2-5 minutes for GitHub Pages to deploy

## Verify Setup

Your documentation should be live at:
```
https://YOUR_USERNAME.github.io/SymbolicDiagonalization.jl/
```

## Automatic Updates

From now on, documentation will automatically update when you:
- Push to `main` branch
- Create a new tag/release

## Testing Before Deployment

Test the build locally before pushing:

```bash
# Simulate CI environment
CI=true julia --project=docs -e '
    using Pkg
    Pkg.develop(PackageSpec(path=pwd()))
    Pkg.instantiate()
'

# Build docs
CI=true julia --project=docs docs/make.jl
```

## Troubleshooting

### "Permission denied (publickey)" Error

**Problem**: SSH key not configured correctly

**Solution**:
1. Verify deploy key is added to GitHub with write access
2. Verify `DOCUMENTER_KEY` secret contains the full private key
3. Check that the key files were generated correctly

### "gh-pages branch not found" Error

**Problem**: First deployment hasn't completed yet

**Solution**:
1. Check Actions tab for deployment status
2. Wait for green checkmark before configuring Pages
3. The branch is created automatically on first successful run

### Documentation Builds but Pages Shows 404

**Problem**: GitHub Pages not configured correctly

**Solution**:
1. Go to Settings → Pages
2. Ensure source is set to `gh-pages` branch
3. Wait 5-10 minutes for propagation
4. Clear browser cache and retry

### Build Fails with "Package not found"

**Problem**: Missing dependencies in `docs/Project.toml`

**Solution**:
1. Add missing packages to `docs/Project.toml`
2. Test locally: `julia --project=docs -e 'using Pkg; Pkg.instantiate()'`

## Security Notes

- **Never commit** `docs/.documenter` (private key) to git
- `.gitignore` is configured to exclude it
- Only the public key goes on GitHub (as a deploy key)
- Only the private key goes in GitHub Secrets

## Need Help?

See the full Contributing Guide in `docs/src/contributing.md` for more details.
