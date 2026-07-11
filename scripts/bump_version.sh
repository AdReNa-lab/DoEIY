#!/usr/bin/env bash
set -e

# Make sure we are at repository root
cd "$(dirname "$0")/.."

VERSION_FILE="version.txt"

if [ ! -f "$VERSION_FILE" ]; then
  echo "Error: $VERSION_FILE not found!" >&2
  exit 1
fi

CURRENT_VERSION=$(cat "$VERSION_FILE" | tr -d '[:space:]')

# Extract major, minor, patch parts
if [[ $CURRENT_VERSION =~ ^([0-9]+)\.([0-9]+)\.([0-9]+)$ ]]; then
  major="${BASH_REMATCH[1]}"
  minor="${BASH_REMATCH[2]}"
  patch="${BASH_REMATCH[3]}"
else
  echo "Error: Current version '$CURRENT_VERSION' in $VERSION_FILE does not match X.Y.Z format!" >&2
  exit 1
fi

BUMP_TYPE=$1

case "$BUMP_TYPE" in
  major)
    major=$((major + 1))
    minor=0
    patch=0
    ;;
  minor)
    minor=$((minor + 1))
    patch=0
    ;;
  patch)
    patch=$((patch + 1))
    ;;
  *)
    echo "Usage: $0 {major|minor|patch}" >&2
    exit 1
    ;;
esac

NEW_VERSION="${major}.${minor}.${patch}"
echo "Bumping version from $CURRENT_VERSION to $NEW_VERSION..."

# Update version.txt
echo "$NEW_VERSION" > "$VERSION_FILE"

# Commit and Tag
git add "$VERSION_FILE"
git commit -m "chore(release): bump version to $NEW_VERSION"
git tag -a "v$NEW_VERSION" -m "Release version $NEW_VERSION"

echo "Successfully bumped version to $NEW_VERSION and created tag v$NEW_VERSION!"
echo "Run 'git push origin main --tags' to push the commit and tag to remote."
