#!/bin/bash

# Organization Name
ORG="cmse822"

# Log files
SUCCESS_LOG="visibility_change_success.log"
FAILURE_LOG="visibility_change_failure.log"

# Initialize log files
echo "Visibility Change Success Log - $(date)" > "$SUCCESS_LOG"
echo "Visibility Change Failure Log - $(date)" > "$FAILURE_LOG"

# Fetch repositories that are both public and archived
REPOS=$(gh repo list "$ORG" --limit 1000 --json nameWithOwner,visibility,isArchived | \
jq -r '.[] | select(.visibility == "PUBLIC" and .isArchived == true) | .nameWithOwner')

# Check if any repositories match the criteria
if [ -z "$REPOS" ]; then
    echo "No public and archived repositories found under the organization '$ORG'."
    exit 0
fi

# Iterate through each repository and change visibility to private
for repo in $REPOS; do
    echo "----------------------------------------"
    echo "Processing repository: $repo"
    
    # Unarchive the repository using GitHub API
    echo "Unarchiving $repo..."
    gh api -X PATCH "/repos/$repo" -F isArchived=false
    if [ $? -eq 0 ]; then
        echo "✅ Successfully unarchived '$repo'." | tee -a "$SUCCESS_LOG"
    else
        echo "❌ Failed to unarchive '$repo'. Skipping visibility change." | tee -a "$FAILURE_LOG"
        echo "----------------------------------------"
        continue
    fi
    
    # Change visibility to private
    echo "Changing visibility of '$repo' to private..."
    gh repo edit "$repo" --visibility private
    if [ $? -eq 0 ]; then
        echo "✅ Successfully changed '$repo' to private." | tee -a "$SUCCESS_LOG"
    else
        echo "❌ Failed to change visibility of '$repo'. Please check your permissions and repository settings." | tee -a "$FAILURE_LOG"
        echo "----------------------------------------"
        continue
    fi
    
    # (Optional) Re-archive the repository if you want to keep it archived
    # Uncomment the following lines if you wish to re-archive after making it private
    echo "Re-archiving '$repo'..."
    gh api -X PATCH "/repos/$repo" -F isArchived=true
    if [ $? -eq 0 ]; then
        echo "✅ Successfully re-archived '$repo'." | tee -a "$SUCCESS_LOG"
    else
        echo "❌ Failed to re-archive '$repo'." | tee -a "$FAILURE_LOG"
    fi

    echo "----------------------------------------"
done

echo "Visibility update process completed."
echo "Check '$SUCCESS_LOG' for successful changes and '$FAILURE_LOG' for any failures."
