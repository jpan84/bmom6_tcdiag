#!/bin/bash -l

set -e

###no trailing slash
ARCHREST="/glade/derecho/scratch/jpan/archive/REPLCASE/rest"
STARTMO="REPLMO"

cd $ARCHREST

if [ "$(pwd)" != "$ARCHREST" ]; then
  echo "Not running in the correct directory. Exiting to prevent rm -rf destroying things."
  exit 1
fi

#### Find all folders that do not match the pattern yy-mm-01-00000
###find . -type d -not -name "*-01-00000" | while read -r folder; do
###  folder_name=$(basename "$folder")
###  if [[ "$folder_name" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}$ ]]; then
###    echo "Removing: $folder"  # Uncomment this line to see what will be removed
###    rm -rf "$folder"        # Uncomment this line to actually remove the folders
###  fi
###done

# Delete all folders, except the 1st of the starting month and every half year thereafter
OTHERMO=$(printf "%02d" $(( ($STARTMO + 6) % 12 )) )
echo $OTHERMO
find . -type d -name "*" | while read -r folder; do
  echo $folder
  folder_name=$(basename "$folder")
  IFS='-' read -r YYYY MM DD SSSSS <<< "$folder_name"
  if [[ "$folder_name" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}$ ]]; then
    if [[ ! "$folder_name" =~ ^[0-9]{4}-$STARTMO-01-00000$ && ! "$folder_name" =~ ^[0-9]{4}-$OTHERMO-01-00000$ ]]; then 
      echo "Removing: $folder"  # Uncomment this line to see what will be removed
      rm -rf "$folder"        # Uncomment this line to actually remove the folders
    fi
  fi
done
