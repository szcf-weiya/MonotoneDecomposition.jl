#!/bin/bash
set -euo pipefail
timestamp=$(date -Iseconds)
competitor=${1:-ss_singe_lambda}
resfolder0=${2:-/tmp}
nlam=${3:-20}
if [[ $competitor == "ss_single_lambda" ]]; then
    nlam=1
fi
nrep=${4:-100}
nfold=${5:-10}
one_se_rule=${6:-false}
# sigmas="[0.1, 0.2, 0.4, 0.5, 1.0, 1.5, 2.0]"
# if [[ $competitor == "ss_fixratio" ]]; then
#     sigmas="0.1:0.1:2.0"
# fi
which julia # in case not found julia, such as run outside the tmux session
for f in "x^2" "x^3" "exp(x)" "sigmoid" "SE_1" "SE_0.1" "Mat12_1" "Mat12_0.1" "Mat32_1" "Mat32_0.1" "RQ_0.1_0.5" "Periodic_0.1_4"; do
#for f in "x^2" "x^3"; do # run locally for debug
    julia ../examples/benchmark.jl $competitor $resfolder0 $nlam $nrep $nfold $one_se_rule $f $timestamp &> nrep$nrep-nfold$nfold-nlam$nlam-$competitor-$f.log &
done
