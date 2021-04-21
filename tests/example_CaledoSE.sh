#!/bin/bash

GONZAG_DATA_DIR="/MEDIA/data/GONZAG"

CLIMPORN_DIR="${HOME}/DEV/climporn/python"

BX="CaledoSE"

F_S="${GONZAG_DATA_DIR}/gonzag_input/SARAL_20170101-20170331.nc"

# CALEDO60 South-Eastern Caledo:
F_M="${GONZAG_DATA_DIR}/gonzag_input/zos_${BX}_2017_JFM_hourly.nc"
V_M="zos"

if [ ! -f ./result.nc ]; then

    CMD="./alongtrack_sat_vs_nemo.py -s ${F_S} -n sla_unfiltered -m ${F_M} -v ${V_M} -l 0"

    echo
    echo "${CMD}"
    ${CMD}
    echo

fi



# Diags:

${CLIMPORN_DIR}/nemo_imshow_2d_field.py ${BX} xnp_msk.nc track 1

${CLIMPORN_DIR}/plot_spectra_SSH_sat_track.py -i result.nc -m ${V_M}_bl -s sla_unfiltered -n 50 \
               -B ${BX} -S "SARAL-AltiKa" -M "CALED60-SouthEast" \
               -a -5 -b 1 -l 13 -L 500
