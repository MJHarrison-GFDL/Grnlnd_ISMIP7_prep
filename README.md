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

e. Run script (python python/create_ISMIP7_inputs.py)

f. Inspect results (Grnlnd_1km.nc, Grnlnd_2km.nc, Grnlnd_4km.nc)
