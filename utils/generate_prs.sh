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

# Automatically set the default repository for GitHub CLI
REPO=$(gh repo view --json nameWithOwner --jq '.nameWithOwner')
if [ -z "$REPO" ]; then
  echo "Error: Could not detect the repository. Are you in the correct directory?"
  exit 1
fi
gh repo set-default "$REPO"

# Ensure the homework directory exists
if [ ! -d "$HOMEWORK_DIR" ]; then
  echo "Error: The homework directory '$HOMEWORK_DIR' does not exist. Please create it first."
  exit 1
fi

# Loop to create notebooks, branches, and PRs
for i in $(seq 1 $NUM_NOTEBOOKS); do
  exercise_name="exercise${i}"
  branch_name="${HOMEWORK_DIR}-exercise-${i}"
  NOTEBOOK_PATH="$HOMEWORK_DIR/${exercise_name}.ipynb"

  # Create and switch to a new branch
  git checkout "$MAIN_BRANCH" || { echo "Error: Failed to checkout $MAIN_BRANCH"; exit 1; }
  # Check if the branch already exists
  if git branch -r | grep "origin/$branch_name" > /dev/null; then
    echo "Branch '$branch_name' already exists. Switching to it."
    git switch "$branch_name"
  else
    # Create the branch
    git switch -C "$branch_name"
  fi

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

  # Create a pull request for the branch
  gh pr create \
    --repo "$REPO" \
    --title "Complete Exercise: $exercise_name" \
    --body "This pull request tracks progress for the **$exercise_name** notebook in Homework **$HOMEWORK_DIR**. Please review and merge when complete." \
    --base "$MAIN_BRANCH" \
    --head "$branch_name" \
    --reviewer "$(gh api repos/$(gh repo view --json nameWithOwner --jq '.nameWithOwner')/collaborators --jq '.[].login' | tr '\n' ',' | sed 's/,$//')"
done

git switch "$MAIN_BRANCH"

echo "All notebooks created, branches pushed, and pull requests opened!"