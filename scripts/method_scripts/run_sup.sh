#!/bin/sh
time -p -f %U_%S bash -c '
python ./methods/superfine/runSuperFine.py -r rmrp "'"$1"'"
'