#!/bin/bash

ffmpeg -i "curl2e-5.1754618_%04d.png" -vf "palettegen=stats_mode=diff" "palette.png"
###ffmpeg -i "curl2e-5.1754618_%04d.png" -vf palettegen -codec:v png palette.png
ffmpeg -framerate 1 -i "./curl2e-5.1754618_%04d.png"  "./output.gif" -i palette.png -lavfi paletteuse
