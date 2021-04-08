#!/bin/bash

GONZAG_DATA_DIR="/MEDIA/data/GONZAG"

CLIMPORN_DIR="${HOME}/DEV/climporn/python"

BX="WestMed"

F_S="${GONZAG_DATA_DIR}/gonzag_input/SENTINEL3A_20170130-20170303.nc"

# eNATL60 WestMed:
F_M="${GONZAG_DATA_DIR}/gonzag_input/sossheig_box_${BX}_eNATL60-BLBT02_20170201-20170228.nc"
V_M="sossheig"

if [ ! -f ./result.nc ]; then

    CMD="./alongtrack_sat_vs_nemo.py -s ${F_S} -n sla_unfiltered -m ${F_M} -v ${V_M} -l 0"
    echo
    echo "${CMD}"
    ${CMD}
    echo

fi

# Diags:

${CLIMPORN_DIR}/nemo_imshow_2d_field.py ${BX} xnp_msk.nc track 1

${CLIMPORN_DIR}/plot_spectra_SSH_sat_track.py -i result.nc -m ${V_M}_bl -s sla_unfiltered -n 70 \
               -B ${BX} -S "Sentinel-3A" -M "eNATL60-WestMed"
