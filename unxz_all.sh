#!/bin/bash

for file in *; do
  if [[ "$file" == *.xz ]]; then
    qcmd -q main -l walltime=00:02:00 -A UPSU0063 unxz -v "${file}"
  fi
done
