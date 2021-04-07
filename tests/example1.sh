#!/bin/bash

GONZAG_DATA_DIR="/MEDIA/data/GONZAG"

CLIMPORN_DIR="${HOME}/DEV/climporn/python"

BX="ORCA1"

#F_S="${GONZAG_DATA_DIR}/gonzag_input/SENTINEL3A_20170801-20170802.nc"
F_S="${GONZAG_DATA_DIR}/gonzag_input/SARAL_20170801-20170801.nc"

F_M="${GONZAG_DATA_DIR}/gonzag_input/ssh_ORCA1_20170101_20171231_grid_T.nc"
V_M="ssh"


if [ ! -f ./result.nc ]; then

    CMD="./alongtrack_sat_vs_nemo.py -s ${F_S} -n sla_unfiltered -m ${F_M} -v ${V_M} -l 0 -D"
    echo
    echo "${CMD}"
    ${CMD}
    echo

fi

# Diags:

${CLIMPORN_DIR}/nemo_imshow_2d_field.py ${BX} xnp_msk.nc track 1

${CLIMPORN_DIR}/plot_spectra_SSH_sat_track.py -i result.nc -m ${V_M}_bl -s sla_unfiltered -B ${BX}
