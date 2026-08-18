#!/bin/bash

module load nco

# Base archive directory
BASE_DIR="/glade/derecho/scratch/jpan/archive/nff_output"

# Path to the mapping file
MAP_FILE="/glade/u/home/jpan/grids/map_ne120np4_to_sx0.66av1_hgrid.stod.250724.nc"

# Loop over all b.e23 directories matching the pattern
for DR in "$BASE_DIR"/b.e23.*; do
    # Check if the path is a valid directory
    if [ -d "$DR" ]; then
        INPUT_DIR="$DR/nff_4mps"
        OUTPUT_DIR="$DR/nff_4mps_sx0.66av1"
        
        # Check if the input directory exists
        if [ -d "$INPUT_DIR" ]; then
            echo "Processing directory: $DR"
            
            # Create the output directory if it doesn't exist
            mkdir -p "$OUTPUT_DIR"
            
            # Loop through all NetCDF files in the input directory
            for FILE in "$INPUT_DIR"/*.h1i.000[6-7]*.nff_4mps; do
                # Check if files exist to handle empty matches safely
                if [ -f "$FILE" ]; then
                    FILENAME=$(basename "$FILE")
                    OUTFILE="$OUTPUT_DIR/$FILENAME"
                    
                    echo "  Regridding: $FILENAME"
                    
                    # Run ncremap for the specific variable TC_R4
                    ncremap -v TC_R4 -m "$MAP_FILE" -i "$FILE" -o "$OUTFILE"
                fi
            done
        else
            echo "Skipping $DR: 'nff_4mps' subfolder not found."
        fi
    fi
done

echo "Regridding complete!"
