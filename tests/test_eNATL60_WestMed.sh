#!/bin/bash

GONZAG_DATA_DIR="/MEDIA/data/GONZAG/gonzag_input"

# eNATL60 WestMed:
F_S="${GONZAG_DATA_DIR}/dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc"

F_M="${GONZAG_DATA_DIR}/sossheig_box_WestMed_eNATL60-BLBT02_20170201-20170228.nc"

./alongtrack_sat_vs_nemo.py -s ${F_S} -n adt_unfiltered \
                            -m ${F_M} -v sossheig \
                            -l 0
