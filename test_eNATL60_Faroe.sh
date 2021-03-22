#!/bin/bash

SOSIE_DIR="${HOME}/DEV/sosie"

F_S="${SOSIE_DIR}/examples/data/dt_global_alg_sla_vxxc_20170402_SARAL-Altika.nc4"

#F_M="${SOSIE_DIR}/examples/data/ssh_ORCA1_20170101_20171231_grid_T.nc4"
#F_L="${SOSIE_DIR}/examples/data/mesh_mask_ORCA1v2_light.nc4"
#ewper="2"

#F_M="${SOSIE_DIR}/examples/data/ssh_ORCA025.nc4"
#F_L="0"
#ewper="2"

#./alongtrack_sat_vs_nemo.py -s ${F_S} -n adt_unfiltered \
#                            -m ${F_M} -v ssh \
#                            -l ${F_L} -p ${ewper}


# eNATL60 Faroe:
#F_S="/MEDIA/data/data/SATELLITE/dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc4"
F_S="/MEDIA/data/data/SATELLITE/SARAL_light_20170107-20170114.nc"
F_M="/MEDIA/data/eNATL60/ZOOMs/Faroe/sossheig_box_Faroe_eNATL60-BLBT02_20170101-20170331_gridT-2D_copy2010.nc4"
F_L="/MEDIA/data/eNATL60/ZOOMs/Faroe/mesh_mask_eNATL60_Faroe.nc"
ewper="-1"

./alongtrack_sat_vs_nemo.py -s ${F_S} -n adt_unfiltered \
                            -m ${F_M} -v sossheig \
                            -l ${F_L} -p ${ewper}
