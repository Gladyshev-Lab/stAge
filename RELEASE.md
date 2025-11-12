# 🚀 stAge Release Creation Instructions

## Automatic Releases via Tags

Releases are created **automatically** when a version tag is created.

### Quick Method (Recommended)

```bash
# Make the script executable (once)
chmod +x release.sh

# Create a release
./release.sh 1.0.0
```

The script automatically:
- ✅ Updates the version in README
- ✅ Creates a release notes template
- ✅ Creates and pushes the tag
- ✅ Triggers the build via GitHub Actions

### Manual Method

```bash
# 1. Ensure all changes are committed
git status

# 2. Create a version tag
git tag -a v1.0.0 -m "Release v1.0.0"

# 3. Push the tag to GitHub
git push origin v1.0.0
```

After this, GitHub Actions will automatically:
- 🏗️ Build the application for Windows, macOS, and Linux
- 📦 Create a release on GitHub
- 📎 Attach the ready binaries

## Version Format

Use **Semantic Versioning**:
- `v1.0.0` - major release (major changes)
- `v1.1.0` - minor release (new features)
- `v1.0.1` - patch release (bug fixes)

## Manual Run via GitHub UI

If you need to rebuild a release without creating a new tag:

1. Go to **Actions** on GitHub
2. Select **Build and Release**
3. Click **Run workflow**
4. Specify the version (e.g., `v1.0.0`)
5. Click **Run workflow**

## Managing Old Releases

By default, only the **last 3 releases** are kept.

To change this:
1. Open `.github/workflows/cleanup-old-releases.yml`
2. Change `keep_latest: 3` to the desired value

To disable auto-deletion:
- Delete the file `.github/workflows/cleanup-old-releases.yml`

## Checking Build Status

```bash
# View the latest releases
gh release list

# View Actions status
gh run list --workflow=build-release.yml
```

## Debugging Issues

### Build Doesn't Start
- ✅ Ensure the tag starts with `v` (e.g., `v1.0.0`)
- ✅ Ensure the tag is pushed: `git push origin v1.0.0`

### Build Fails
- 🔍 Check logs in GitHub Actions
- 🔍 Ensure all dependencies are in `requirements.txt`
- 🔍 Check permissions (Settings → Actions → General → Workflow permissions)

### Release Not Created
- ✅ Ensure write permissions are enabled: Settings → Actions → General → Workflow permissions → Read and write permissions

## Local Build for Testing

Before creating a release, you can test the build locally:

```bash
# macOS
pyinstaller --name="stAge" --windowed --onefile main.py
open dist/stAge.app

# Test run
./dist/stAge.app/Contents/MacOS/stAge
```

## Pre-Release Checklist

- [ ] All changes are committed
- [ ] Tests passed (if any)
- [ ] README.md updated
- [ ] Version updated
- [ ] CHANGELOG filled
- [ ] Local build works

## Example Commands

```bash
# Create a patch release
./release.sh 1.0.1

# Create a minor release
./release.sh 1.1.0

# Create a major release
./release.sh 2.0.0

# Delete a tag (if mistaken)
git tag -d v1.0.0
git push origin :refs/tags/v1.0.0

# Recreate a tag
git tag -a v1.0.0 -m "Release v1.0.0"
git push origin v1.0.0
```

## Release Structure

Each release includes:
- 📦 **Windows**: `stAge-Windows-x64.zip` (~50-100 MB)
- 📦 **macOS**: `stAge-macOS-ARM64.dmg` (~80-120 MB)
- 📦 **Linux**: `stAge-Linux-x64.tar.gz` (~70-110 MB)
- 📝 **Release Notes**: description of changes
- 🔗 **Source Code**: source code in zip/tar.gz

## GitHub CLI (Optional)

For advanced management, install GitHub CLI:

```bash
# macOS
brew install gh

# Authentication
gh auth login

# Create a release from scratch
gh release create v1.0.0 \
  --title "stAge v1.0.0" \
  --notes "Release notes here"

# View releases
gh release list

# Delete a release
gh release delete v1.0.0
```