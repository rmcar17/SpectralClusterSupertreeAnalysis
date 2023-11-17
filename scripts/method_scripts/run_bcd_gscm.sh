#!/bin/sh
out_path=${1//[\/.]/_}_tmp_bcd_gscm.nwk
time -p -f %U_%S bash -c '
java -jar ./methods/bcd/BCDSupertrees.jar -L OFF -s true -w BRANCH_LENGTH -f NEWICK -d NEWICK -o '"$out_path"' -t 1 -B "'"$1"'" &> /dev/null
wait # weighting (-w) can either be UNIT_WEIGHT, TREE_WEIGHT, BRANCH_LENGTH or mroe
cat '"$out_path"'
wait
rm '"$out_path"'
wait
'