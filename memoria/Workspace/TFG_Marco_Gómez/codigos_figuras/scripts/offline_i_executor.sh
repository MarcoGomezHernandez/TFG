#!/bin/bash

BASE_DIR="data"
FILE_PREFIX="best_sample_syn_params_i"
FILE_BASE_PATH="${BASE_DIR}/${FILE_PREFIX}"
HISTORY_PREFIX="bo_history_offline_syn_i"
HISTORY_BASE_PATH="${BASE_DIR}/${HISTORY_PREFIX}"

CMD_BASE="../codigos/offline/build/bo ../codigos/offline/PDtracesV.csv 0 0.008 800 320 560 560 50 350"
CMD_END="-4.0 4.0 1"

PHASE_ARGS="1 0 1 2 0.3 -1.4 -0.4"
ANTIPHASE_ARGS="0 0 1 2 0.3 0.2 0.8"

run_optimization() {
  local display_name=$1
  local mode_name=$2
  local mode_args=$3
  
  echo "$display_name"
  local b_s=""
  local b_f=""
  
  for i in {0..19}; do
    echo "Muestra $i de 20"
    local f="${FILE_BASE_PATH}_${mode_name}_${i}.yaml"
    local h="${HISTORY_BASE_PATH}_${mode_name}.jsonl"
    
    local o=$($CMD_BASE $mode_args $CMD_END "$f" "$h")
    local s=$(awk '/Best:/{print $2}' <<< "$o")
    
    if [[ $i -eq 1 ]] || awk "BEGIN{exit !($s > $b_s)}"; then 
      b_s=$s
      b_f=$f
    fi
  done
}

run_optimization "Fase" "phase" "$PHASE_ARGS"
run_optimization "Antifase" "antiphase" "$ANTIPHASE_ARGS"