#!/bin/bash

BASE_DIR="data"

CMD_BASE="../codigos/offline/build/bo ../codigos/offline/PDtracesV.csv 0 0.008 800 320 560 560 50 350"
CMD_END="-4.0 4.0 1"

PHASE_ARGS="1 0 1 2 0.3 -1.4 -0.4"
ANTIPHASE_ARGS="0 0 1 2 0.3 0.2 0.8"

run_optimization() {
  local display_name=$1
  local mode_name=$2
  local mode_args=$3
  
  echo "$display_name"
  
  for i in {0..19}; do
    echo "Muestra $i de 20"
    local f="${BASE_DIR}/sample_syn_params_i_${mode_name}_${i}.yaml"
    local h="${BASE_DIR}/bo_history_offline_syn_i_${mode_name}.jsonl"

    $CMD_BASE $mode_args $CMD_END "$f" "$h" >/dev/null 2>&1
  done
}

run_optimization "Fase" "phase" "$PHASE_ARGS"
run_optimization "Antifase" "antiphase" "$ANTIPHASE_ARGS"