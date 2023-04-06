#!/bin/bash

source pretty_latlon.sh


output_dir=output
beg_year=2016
end_year=2018


# format: lat_m lat_M lon_m lon_M
spatial_rngs=(
    37 38 -123 -122
    31 43 -130 -120
)


nparms=4

for i in $( seq 1 $(( ${#spatial_rngs[@]} / $nparms )) ); do

    lat_min=${spatial_rngs[$(( ( i - 1 ) * $nparms + 0 ))]}
    lat_max=${spatial_rngs[$(( ( i - 1 ) * $nparms + 1 ))]}
    lon_min=${spatial_rngs[$(( ( i - 1 ) * $nparms + 2 ))]}
    lon_max=${spatial_rngs[$(( ( i - 1 ) * $nparms + 3 ))]}

    time_str=$( printf "%04d-%04d" $beg_year $end_year )
    spatial_str=$( printf "%s-%s_%s-%s" $( pretty_lat $lat_min ) $( pretty_lat $lat_max ) $( pretty_lon $lon_min ) $( pretty_lon $lon_max ) )


    output_file="$output_dir/ARinabox_${time_str}_${spatial_str}.nc"

    echo "time_str    : $time_str"    
    echo "spatial_str : $spatial_str"
    echo "output_file : $output_file"

    mkdir -p $output_dir

    eval "python3 AR_presense_in_a_box.py \\
        --year-rng $beg_year $end_year \\
        --lat $lat_min $lat_max \\
        --lon $lon_min $lon_max \\
        --output-file $output_file
    " &

    wait

done
