#!/bin/bash


# eNATL60 Faroe:
#F_S="/MEDIA/data/data/SATELLITE/dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc4"
F_S="/MEDIA/data/data/SATELLITE/SARAL_light_20170107-20170114.nc"

F_M="/MEDIA/data/eNATL60/ZOOMs/Faroe/sossheig_box_Faroe_eNATL60-BLBT02_20170101-20170331.nc"
F_L="/MEDIA/data/eNATL60/ZOOMs/Faroe/mesh_mask_eNATL60_Faroe.nc"
ewper="-1"

./alongtrack_sat_vs_nemo.py -s ${F_S} -n adt_unfiltered \
                            -m ${F_M} -v sossheig \
                            -l ${F_L} -p ${ewper}
