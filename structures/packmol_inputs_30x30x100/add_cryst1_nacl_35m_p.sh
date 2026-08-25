#!/bin/bash
# insert the correct CRYST1 record into every packmol output
CRYST="CRYST1   30.000   30.000  100.000  90.00  90.00  90.00 P 1           1"
for f in nacl_35m_pr*.pdb; do
  grep -q "^CRYST1" "$f" && sed -i "/^CRYST1/d" "$f"
  sed -i "1i $CRYST" "$f"
  echo "$f: $(grep -c HOH "$f") HOH  $(head -1 "$f")"
done
