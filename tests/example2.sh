#!/bin/bash

GONZAG_DATA_DIR="/MEDIA/data/GONZAG"

CLIM_PORN_DIR="${HOME}/DEV/climporn/python"

BX="SouthWestPac_G12"

F_S="${GONZAG_DATA_DIR}/gonzag_input/SARAL_20190301-20190331.nc"

F_M="${GONZAG_DATA_DIR}/gonzag_input/zos_GLORYS12V1_SouthWestPac_20190301-20190331_hourly.nc"
V_M="zos"


if [ ! -f ./result.nc ]; then

    CMD="./alongtrack_sat_vs_nemo.py -s ${F_S} -n sla_unfiltered -m ${F_M} -v ${V_M} -l 0"
    echo
    echo "${CMD}"
    ${CMD}
    echo

fi

# Diags:

${CLIM_PORN_DIR}/nemo_imshow_2d_field.py ${BX} xnp_msk.nc track 1

${CLIM_PORN_DIR}/plot_spectra_SSH_sat_track.py -i result.nc -m ${V_M}_bl -s sla_unfiltered -B ${BX}
