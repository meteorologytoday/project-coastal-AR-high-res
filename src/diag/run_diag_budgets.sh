#!/bin/bash


output_dir="output_diag_budgets"
data_dir="/data/SO2/SWOT/MARA/RUN4_LY/DIAGS_DLY"
grid_dir="/data/SO2/SWOT/GRID/BIN"
beg_date="2015-08-15"
end_date="2018-08-22"
mitgcm_beg_date="2015-01-01"
mitgcm_deltaT=150.0
nproc=16

python3 diag_budgets.py \
    --data-dir $data_dir \
    --grid-dir $grid_dir \
    --mitgcm-beg-date $mitgcm_beg_date \
    --mitgcm-deltaT $mitgcm_deltaT \
    --beg-date $beg_date \
    --end-date $end_date \
    --nproc $nproc \
    --output-dir $output_dir




