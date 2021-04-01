#!/bin/bash

GONZAG_DATA_DIR="/MEDIA/data/GONZAG/gonzag_input"

# eNATL60 Faroe:
F_S="${GONZAG_DATA_DIR}/dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc"

F_M="${GONZAG_DATA_DIR}/sossheig_box_Faroe_eNATL60-BLBT02_20170101-20170331.nc"

F_L=${F_M}
V_L="tmask"

#ewper="-1"

./alongtrack_sat_vs_nemo.py -s ${F_S} -n adt_unfiltered \
                            -m ${F_M} -v sossheig \
                            -l ${F_L} -k ${V_L} -D
