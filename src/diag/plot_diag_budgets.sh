#!/bin/bash

avg_days=3
diag_data_dir="output_diag_budgets"
figure_dir="figures/diag_budgets_${avg_days}"
beg_date="2016-10-10"
#end_date="2017-03-31"
end_date="2016-10-10"

python3 plot_diag_budgets.py \
    --date-rng $beg_date $end_date \
    --input-dir $diag_data_dir \
    --output-dir $figure_dir \
    --avg-days $avg_days  \
    --varnames heat.dMLTdt heat.G_sw_nsw heat.G_adv_fwf heat.G_vdiff_ent \
               salt.dMLSdt salt.G_sfc    salt.G_adv_fwf salt.G_vdiff_ent \
               heat.MLD    heat.dMLDdt \
    --ncol 4 \
    --nproc 1 --overwrite

if [ ] ; then
python3 plot_diag_budgets.py \
    --date-rng $beg_date $end_date \
    --input-dir $diag_data_dir \
    --output-dir $figure_dir \
    --avg-days $avg_days  \
    --varnames heat.dMLTdt heat.G_sw  heat.G_nsw  heat.G_adv  heat.G_fwf  heat.G_vdiff  heat.G_ent \
               salt.dMLSdt BLANK      salt.G_sfc  salt.G_adv  salt.G_fwf  salt.G_vdiff  salt.G_ent \
               heat.MLD    heat.dMLDdt \
    --ncol 7 \
    --nproc 1
fi

