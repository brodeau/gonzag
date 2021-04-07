#!/bin/bash

GONZAG_DATA_DIR="/MEDIA/data/GONZAG"

CLIM_PORN_DIR="${HOME}/DEV/climporn/python"

BX="Faroe"

F_S="${GONZAG_DATA_DIR}/gonzag_input/SARAL_20170101-20170331.nc"

# eNATL60 Faroe:
F_M="${GONZAG_DATA_DIR}/gonzag_input/sossheig_box_${BX}_eNATL60-BLBT02_20170101-20170331.nc"

if [ ! -f ./result.nc ]; then

    CMD="./alongtrack_sat_vs_nemo.py -s ${F_S} -n sla_unfiltered -m ${F_M} -v sossheig \
                                     -l ${F_M} -k tmask -D"

    echo
    echo "${CMD}"
    ${CMD}
    echo

fi


# Diags:

${CLIM_PORN_DIR}/nemo_imshow_2d_field.py ${BX} xnp_msk.nc track 1

${CLIM_PORN_DIR}/plot_spectra_SSH_sat_track.py -i result.nc -m sossheig_bl -s sla_unfiltered -B ${BX}
