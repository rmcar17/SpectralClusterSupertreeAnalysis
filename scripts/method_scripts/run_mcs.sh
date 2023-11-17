#!/bin/sh
time -p -f %U_%S bash -c '
python ./methods/mcs/run_mcs.py "'"$1"'"
'