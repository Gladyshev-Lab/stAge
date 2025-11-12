#!/bin/bash

# Script for creating a new stAge release
# Usage: ./release.sh 1.0.0

set -e

if [ -z "$1" ]; then
    echo "❌ Error: specify the version"
    echo "Usage: ./release.sh 1.0.0"
    exit 1
fi

VERSION="$1"
TAG="v${VERSION}"

echo "🚀 Creating stAge release ${TAG}"
echo ""

# Check that we are in a git repository
if ! git rev-parse --git-dir > /dev/null 2>&1; then
    echo "❌ Error: git repository not found"
    exit 1
fi

# Check that there are no uncommitted changes
if ! git diff-index --quiet HEAD --; then
    echo "⚠️  Warning: there are uncommitted changes"
    read -p "Continue? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Update version in README.md (if badge exists)
if [ -f "README.md" ]; then
    echo "📝 Updating version in README.md..."
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS
        sed -i '' "s/v[0-9]\+\.[0-9]\+\.[0-9]\+/${TAG}/g" README.md
    else
        # Linux
        sed -i "s/v[0-9]\+\.[0-9]\+\.[0-9]\+/${TAG}/g" README.md
    fi
fi

# Create CHANGELOG for this release
echo "📋 Generating CHANGELOG..."
cat > RELEASE_NOTES.md << EOF
# Release ${TAG}

## 🎉 What's New

- [Add new features here]

## 🐛 Bug Fixes

- [Add fixed bugs here]

## 🔧 Improvements

- [Add improvements here]

---

**Full Changelog**: https://github.com/$(git config --get remote.origin.url | sed 's/.*github.com[:/]\(.*\)\.git/\1/')/compare/v${VERSION}

EOF

# Open in editor for editing
if command -v code &> /dev/null; then
    code RELEASE_NOTES.md
elif command -v nano &> /dev/null; then
    nano RELEASE_NOTES.md
else
    vi RELEASE_NOTES.md
fi

echo ""
read -p "📝 Edit RELEASE_NOTES.md and press Enter to continue..."

# Commit changes
echo "💾 Committing changes..."
git add README.md RELEASE_NOTES.md
git commit -m "Release ${TAG}" || true

# Create tag
echo "🏷️  Creating tag ${TAG}..."
if git rev-parse "$TAG" >/dev/null 2>&1; then
    echo "⚠️  Tag ${TAG} already exists!"
    read -p "Delete existing tag? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        git tag -d "$TAG"
        git push origin :refs/tags/"$TAG" || true
    else
        exit 1
    fi
fi

git tag -a "$TAG" -m "Release ${TAG}"

# Push to GitHub
echo "📤 Pushing to GitHub..."
git push origin main
git push origin "$TAG"

echo ""
echo "✅ Done!"
echo ""
echo "🔗 GitHub Actions will start building automatically"
echo "📦 Release will appear here: https://github.com/$(git config --get remote.origin.url | sed 's/.*github.com[:/]\(.*\)\.git/\1/')/releases"
echo ""
echo "⏳ Wait ~10-15 minutes for completion of builds for all platforms"