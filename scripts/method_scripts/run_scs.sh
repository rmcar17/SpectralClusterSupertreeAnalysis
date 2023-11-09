#!/bin/sh
export PYTHONHASHSEED=0
time python ./methods/scs/run_scs.py $1 $2 "$3"