#!/bin/bash

set -e

# Default remote and paths
REMOTE=origin
ARCHIVE_BRANCH=archive
PUBLISH_BRANCH=gh-pages
ARCHIVE_DIR=archive

# Make sure we are not in a dirty state
git status --porcelain | grep . && { echo "Please commit or stash changes before running this script."; exit 1; }

# Fetch latest branches
git fetch $REMOTE $ARCHIVE_BRANCH
git fetch $REMOTE $PUBLISH_BRANCH

# Switch to gh-pages branch
git checkout $PUBLISH_BRANCH
git pull $REMOTE $PUBLISH_BRANCH

# Copy archive folder from archive branch
git checkout $ARCHIVE_BRANCH -- $ARCHIVE_DIR

# Stage, commit, and push if there are changes
if git status --porcelain | grep $ARCHIVE_DIR; then
  git add $ARCHIVE_DIR
  git commit -m "Update $ARCHIVE_DIR from $ARCHIVE_BRANCH"
  git push $REMOTE $PUBLISH_BRANCH
  echo "Archive folder updated on $PUBLISH_BRANCH."
else
  echo "No changes in $ARCHIVE_DIR to commit."
fi

# Switch back to previous branch (optional)
git checkout -
