#!/bin/bash

GONZAG_DATA_DIR="/MEDIA/data/GONZAG"

CLIM_PORN_DIR="${HOME}/DEV/climporn/python"

BX="WestMed"

F_S="${GONZAG_DATA_DIR}/gonzag_input/SENTINEL3A_20170130-20170303.nc"

# eNATL60 WestMed:
F_M="${GONZAG_DATA_DIR}/gonzag_input/sossheig_box_${BX}_eNATL60-BLBT02_20170201-20170228.nc"

./alongtrack_sat_vs_nemo.py -s ${F_S} -n sla_unfiltered \
                            -m ${F_M} -v sossheig \
                            -l 0

# Diags:

${CLIM_PORN_DIR}/nemo_imshow_2d_field.py ${BX} xnp_msk.nc track 1

${CLIM_PORN_DIR}/plot_spectra_SSH_sat_track.py -i result.nc -m sossheig_bl -s sla_unfiltered -B ${BX}
