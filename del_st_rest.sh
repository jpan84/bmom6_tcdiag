#!/bin/bash

ARCHREST="/glade/derecho/scratch/archive/jpan/CASENAME/rest"

cd $ARCHREST

# Find all folders that do not match the pattern yy-mm-01-00000
find . -type d -not -name "*-01-00000" | while read -r folder; do
  echo "Removing: $folder"  # Uncomment this line to see what will be removed
  rm -rf "$folder"        # Uncomment this line to actually remove the folders
done
