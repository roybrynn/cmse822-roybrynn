#!/bin/bash

# Check for required inputs
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <target_directory> <number_of_exercises>"
  echo "Example: $0 homework1 10"
  exit 1
fi

# Command-line inputs
HOMEWORK_DIR=$1
NUM_NOTEBOOKS=$2

# Configuration
MAIN_BRANCH="main" # The branch to base all new branches on

# Ensure GitHub CLI (gh) is installed
if ! command -v gh &> /dev/null; then
  echo "Error: GitHub CLI (gh) is not installed. Please install it before running this script."
  exit 1
fi

# Detect the current repository based on the local Git remote URL
REPO=$(git config --get remote.origin.url | sed -E 's#(git@|https://)([^:/]+)[:/]([^/]+)/([^/.]+)(\.git)?#\3/\4#')
if [ -z "$REPO" ]; then
  echo "Error: Could not detect the repository. Are you in the correct directory?"
  exit 1
fi

# Set the detected repository as the default for the GitHub CLI
gh repo set-default "$REPO"

# Ensure the homework directory exists
mkdir -p "$HOMEWORK_DIR"

# Loop to create notebooks, branches, and PRs
for i in $(seq 1 $NUM_NOTEBOOKS); do
  exercise_name="exercise${i}"
  branch_name="${HOMEWORK_DIR}-exercise-${i}"
  NOTEBOOK_PATH="$HOMEWORK_DIR/${exercise_name}.ipynb"

  # Create and switch to a new branch
  git checkout "$MAIN_BRANCH" || { echo "Error: Failed to checkout $MAIN_BRANCH"; exit 1; }
  git switch -c "$branch_name"

  # Create the notebook JSON structure
  cat <<EOF > "$NOTEBOOK_PATH"
{
  "cells": [],
  "metadata": {},
  "nbformat": 4,
  "nbformat_minor": 4
}
EOF
  echo "Generated $NOTEBOOK_PATH"

  # Commit and push the changes
  git add "$NOTEBOOK_PATH"
  git commit -m "Add $exercise_name notebook"
  git push origin "$branch_name"

  # Reviewers
    echo "$(gh api repos/$REPO/collaborators --jq '.[].login' | tr '\n' ',' | sed 's/,$//')"

  # Create a pull request for the branch
  gh pr create \
    --title "Complete Exercise: $exercise_name" \
    --body "This pull request tracks progress for the **$exercise_name** notebook in Homework **$HOMEWORK_DIR**. Please review and merge when complete." \
    --base "$MAIN_BRANCH" \
    --head "$branch_name" \
    --reviewer "$(gh api repos/$REPO/collaborators --jq '.[].login' | tr '\n' ',' | sed 's/,$//')"
done

git checkout "$MAIN_BRANCH"

echo "All notebooks created, branches pushed, and pull requests opened!"
