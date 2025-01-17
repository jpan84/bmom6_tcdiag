#!/bin/bash

set -e

ARCHREST="/glade/derecho/scratch/jpan/archive/CASENAME/rest"

cd $ARCHREST

if [ "$(pwd)" != "$ARCHREST" ]; then
  echo "Not running in the correct directory. Exiting to prevent rm -rf destroying things."
  exit 1
fi

# Find all folders that do not match the pattern yy-mm-01-00000
find . -type d -not -name "*-01-00000" | while read -r folder; do
  folder_name=$(basename "$folder")
  if [[ "$folder_name" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}$ ]]; then
    echo "Removing: $folder"  # Uncomment this line to see what will be removed
    #rm -rf "$folder"        # Uncomment this line to actually remove the folders
  fi
done