This repository is used to generate needed input files for MOM6 Greenland ISMIP7 experiments

Reference:
==========

https://www.ismip.org/home

Requirements:
=============

Python (version 3)



Instructions:
=============

a. Create a anaconda environement:
   	  conda env create -f conda_environment.yml

b. Activate conda environment:
   	  conda activate prep_ISMIP7

c. Download ISMIP7 datasets (e.g. /path_to_ISMIP_data)

d. Link to current directory  ( ln -s /path_to_ISMIP_data INPUT)

e.  ln -s /path_to_rheology INPUT_NS

f. Run scripts

   - python python/Grnld_bed_IS.py
   - python python/Grnlnd_forcing.py
   - python python/Grnlnd_rheology.py



g. Replace files with updated hmask (this step will go away)

   - mv Grnld_hmask_1km.nc Grnld_1km.nc
   - mv Grnld_hmask_2km.nc Grnld_2km.nc
   - mv Grnld_hmask_4km.nc Grnld_4km.nc
