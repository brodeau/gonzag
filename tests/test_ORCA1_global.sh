#!/bin/bash

GONZAG_DATA_DIR="/MEDIA/data/gonzag_input"

F_S="${GONZAG_DATA_DIR}/dt_global_alg_sla_vxxc_20170402_SARAL-Altika.nc"
#F_S="${GONZAG_DATA_DIR}/saral_short_20170402.nc"

F_M="${GONZAG_DATA_DIR}/ssh_ORCA1_20170101_20171231_grid_T.nc"

F_L="0"

#ewper="2"

./alongtrack_sat_vs_nemo.py -s ${F_S} -n adt_unfiltered \
                            -m ${F_M} -v ssh \
                            -l ${F_L} -D


