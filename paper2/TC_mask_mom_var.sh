#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e
module load cdo

# Base output directory where regridded masks reside
MASK_BASE="/glade/derecho/scratch/jpan/archive/nff_output"

# Base directory for the ocean history files
OCN_BASE="/glade/campaign/univ/upsu0032/jpan_aquaptc"

for DR_PATH in "$MASK_BASE"/b.e23.*; do
    if [ -d "$DR_PATH" ]; then
        DR=$(basename "$DR_PATH")
        MASK_DIR="$DR_PATH/nff_4mps_sx0.66av1"
        OCN_DIR="$OCN_BASE/$DR/ocn/hist"
        
        # Verify both required directories exist
        if [ -d "$MASK_DIR" ] && [ -d "$OCN_DIR" ]; then
            echo "Processing experiment: $DR"
            
            # Temporary directory for intermediate files
            TMP_DIR="$DR_PATH/tmp_cdo_processing"
            mkdir -p "$TMP_DIR"
            
            MASK_CAT="$TMP_DIR/cat_mask.nc"
            OCN_CAT="$TMP_DIR/cat_ocn.nc"
            OCN_SUBSET="$TMP_DIR/subset_ocn.nc"
            FINAL_OUTPUT="$DR_PATH/hflso_masked_4mps.nc"
            
            # Define output file paths for zonal means
            OCN_SUBSET_ZONAL="$DR_PATH/hflso_subset_zonal_mean.nc"
            FINAL_ZONAL="$DR_PATH/hflso_masked_4mps_zonal_mean.nc"
            
            echo "  1) Concatenating mask files..."
            cdo cat "$MASK_DIR"/*.nff_4mps "$MASK_CAT"
            
            echo "  2) Concatenating ocean files and matching timesteps..."
            cdo cat "$OCN_DIR"/*.sfc.*.nc "$OCN_CAT"
            
            # Select hflso and match the timesteps present in the concatenated mask
            cdo selname,hflso -selyearmondayhour,$(cdo showtimestamp "$MASK_CAT" | tr ' ' '\n' | grep -v '^$' | tr '\n' ',') "$OCN_CAT" "$OCN_SUBSET"
            
            echo "  3) Multiplying hflso by mask..."
            # Uses element-wise multiplication (mul) for matching spatial/temporal dimensions
            cdo mul "$OCN_SUBSET" "$MASK_CAT" "$FINAL_OUTPUT"
            
            echo "  4) Computing zonal means..."
            # Zonal mean of the subsetted ocean output
            cdo zonmean "$OCN_SUBSET" "$OCN_SUBSET_ZONAL"
            
            # Zonal mean of the final masked output
            cdo zonmean "$FINAL_OUTPUT" "$FINAL_ZONAL"
            
            # Clean up intermediate directory
            rm -rf "$TMP_DIR"
            echo "  Completed -> Outputs saved to:"
            echo "    - 2D Masked Output:  $FINAL_OUTPUT"
            echo "    - Ocn Zonal Mean:    $OCN_SUBSET_ZONAL"
            echo "    - Masked Zonal Mean: $FINAL_ZONAL"
        else
            echo "Skipping $DR: Required mask or ocean directory not found."
        fi
    fi
done

echo "All calculations complete!"
