#!/bin/bash

avg_days=3
diag_data_dir="output_diag_budgets"
figure_dir="figures/diag_budgets_${avg_days}"
beg_date="2016-10-01"
end_date="2017-03-31"

python3 plot_diag_budgets.py \
    --date-rng $beg_date $end_date \
    --input-dir $diag_data_dir \
    --output-dir $figure_dir \
    --avg-days $avg_days  \
    --no-display


