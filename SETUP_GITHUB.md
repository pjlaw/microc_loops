# Setup Instructions for GitHub

## Option 1: Create Repository on GitHub (Recommended)

If the `plaw/microc_loops` repository doesn't exist yet:

### Step 1: Create the Repository

1. Go to https://github.com/new
2. Set repository name: `microc_loops`
3. Owner: `plaw`
4. Description: "Hi-C loop detection pipeline using raichu and pyHICCUPS"
5. Choose Public or Private
6. **Do NOT initialize with README, .gitignore, or license** (we already have these)
7. Click "Create repository"

### Step 2: Push Your Code

From the `microc_loops` directory on your local machine:

```bash
cd microc_loops

# Make sure you're on the right branch
git branch

# Add the remote (if not already added)
git remote add origin https://github.com/plaw/microc_loops.git

# Push the branch
git push -u origin seqera-ai/20260402-082858-initial-pipeline

# Or push to main/master
git checkout -b main  # or master
git push -u origin main
```

### Step 3: Create Pull Request (Optional)

If you pushed to a feature branch:

1. Go to https://github.com/plaw/microc_loops
2. Click "Compare & pull request"
3. Review changes
4. Click "Create pull request"

---

## Option 2: Repository Already Exists

If `plaw/microc_loops` already exists but is empty or you want to add these files:

### Clone and Add Files

```bash
# Clone the existing repository
git clone https://github.com/plaw/microc_loops.git
cd microc_loops

# Copy the pipeline files (from wherever you extracted them)
cp -r /path/to/pipeline/files/* .

# Create a new branch
git checkout -b seqera-ai/20260402-082858-initial-pipeline

# Add and commit
git add .
git commit -m "Initial commit: Hi-C loop detection pipeline with raichu and pyHICCUPS"

# Push
git push -u origin seqera-ai/20260402-082858-initial-pipeline
```

---

## Option 3: Using the Archive

Download the `microc_loops.tar.gz` archive and extract it:

```bash
# Extract
tar -xzf microc_loops.tar.gz
cd microc_loops

# The repository is already initialized with all commits
# Just add the remote and push
git remote add origin https://github.com/plaw/microc_loops.git
git push -u origin seqera-ai/20260402-082858-initial-pipeline
```

---

## Verify the Push

After pushing, verify on GitHub:

1. Visit https://github.com/plaw/microc_loops
2. Check that the branch appears in the branch dropdown
3. Verify all files are present:
   - main.nf
   - nextflow.config
   - README.md
   - .gitignore

---

## Files Included

✅ **main.nf** - Complete Nextflow pipeline (312 lines)
  - Raichu normalization
  - pyHICCUPS loop calling at multiple resolutions
  - Loop combination and format conversion
  - Statistics generation

✅ **nextflow.config** - Full configuration (146 lines)
  - Multiple execution profiles (docker, singularity, conda, slurm, aws, gcp)
  - Resource specifications
  - Reporting configuration

✅ **README.md** - Comprehensive documentation (228 lines)
  - Quick start guide
  - Parameter descriptions
  - Output structure
  - Troubleshooting

✅ **.gitignore** - Proper exclusions
  - Nextflow work directories
  - Data files
  - Reports and logs

---

## Troubleshooting

### Authentication Issues

If you get authentication errors:

```bash
# Use personal access token
git remote set-url origin https://YOUR_TOKEN@github.com/plaw/microc_loops.git
git push -u origin seqera-ai/20260402-082858-initial-pipeline
```

### Repository Not Found

Make sure:
1. The repository exists at https://github.com/plaw/microc_loops
2. You have write access to it
3. The repository name is spelled correctly

### Permission Denied

Check that:
1. You're logged in to GitHub
2. You have permissions for the `plaw` account/organization
3. Your SSH key or token is configured correctly

---

## Need Help?

If you encounter issues:
1. Check the GitHub repository settings
2. Verify your authentication (SSH key or personal access token)
3. Ensure you have write permissions to the repository
