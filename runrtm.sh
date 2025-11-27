#!/bin/bash

#SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
#SCRIPT_DIR=$(dirname "$(realpath "$0")")
#PAR_FILE="$SCRIPT_DIR/par/parameters.par"

#if [ ! -f "$PAR_FILE" ]; then
#    echo "Error: Cannot find parameter file: $PAR_FILE"
#    exit 1
#fi

# forward model
#./bin/fd "par=$PAR_FILE"
# RTM
#./bin/rtm "par=$PAR_FILE"

export CUDA_VISIBLE_DEVICES=1

# forward model
#(cd "$SCRIPT_DIR/par" && ../bin/fd "par=./parameters.par")
# RTM
#(cd "$SCRIPT_DIR/par" && ../bin/rtm "par=./parameters.par")



nohup ./bin/upc_run par=./rtm.par > ./upc_runrtm.log 2>&1 &
