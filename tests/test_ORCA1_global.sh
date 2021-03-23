#!/bin/bash

SOSIE_DIR="${HOME}/DEV/sosie"

F_S="${SOSIE_DIR}/examples/data/dt_global_alg_sla_vxxc_20170402_SARAL-Altika.nc4"

F_M="${SOSIE_DIR}/examples/data/ssh_ORCA1_20170101_20171231_grid_T.nc4"
#F_L="${SOSIE_DIR}/examples/data/mesh_mask_ORCA1v2_light.nc4"
F_L="0"
ewper="2"

#F_M="${SOSIE_DIR}/examples/data/ssh_ORCA025.nc4"
#F_L="0"
#ewper="2"

./alongtrack_sat_vs_nemo.py -s ${F_S} -n adt_unfiltered \
                            -m ${F_M} -v ssh \
                            -l ${F_L} -p ${ewper}
