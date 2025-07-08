###Concatenate the logs in CESM archive
#!/bin/bash
comps=("atm" "cesm" "diags" "drv" "ice" "lnd" "med" "ocn")

for cmp in "${comps[@]}"; do
  echo "catting ${cmp} logs"
  find . -name "${cmp}.log.*.gz" -print0 | sort -z | xargs -0 cat > ${cmp}.log.cat.gz
done
