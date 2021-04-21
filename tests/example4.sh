#!/bin/bash

GONZAG_DATA_DIR="/MEDIA/data/GONZAG"

CLIMPORN_DIR="${HOME}/DEV/climporn/python"

BX="Faroe"

F_S="${GONZAG_DATA_DIR}/gonzag_input/SARAL_20170101-20170331.nc"

# eNATL60 Faroe:
F_M="${GONZAG_DATA_DIR}/gonzag_input/sossheig_box_${BX}_eNATL60-BLBT02_20170101-20170331.nc"
V_M="sossheig"

if [ ! -f ./result.nc ]; then

    CMD="./alongtrack_sat_vs_nemo.py -s ${F_S} -n sla_unfiltered -m ${F_M} -v ${V_M} \
                                     -l ${F_M} -k tmask -D"

    echo
    echo "${CMD}"
    ${CMD}
    echo

fi

# Diags:

${CLIMPORN_DIR}/nemo_imshow_2d_field.py ${BX} xnp_msk.nc track 1

${CLIMPORN_DIR}/plot_spectra_SSH_sat_track.py -i result.nc -m ${V_M}_bl -s sla_unfiltered -n 60 -B ${BX} -S "Saral-Altika" -M "eNATL60-Faroe" \
                                              -a -3 -b 1 -l 13 -L 500
