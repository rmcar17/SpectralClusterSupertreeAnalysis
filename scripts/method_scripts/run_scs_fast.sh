#!/bin/sh
export PYTHONHASHSEED=0
time -p -f %U_%S bash -c '
python ./methods/scs/run_scs.py -c '"$1"' '"$2"' "'"$3"'"
'