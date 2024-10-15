Code to reconstruct vegetation from a subset of pollen records following the LegacyVegetation
dataset. In this repository parts of the analysis can be replicated using a subset of three records.

The script main_complete.R includes all of the analyses and takes approximately 5 minutes to run all
scripts in the necessary order. It is recommended to run this first as it also takes care of packages
that need to be installed.

The following scripts are included in this analysis:
	- main_REVEALS.R
		reconstructs vegetation cover from three pollen records using the REVEALS model and the
		original RPP estimates
		output: output/PANGAEA/Europe_original_RPP.csv
	- main_reconstruct_forest.R
		reconstruct forests cover from Pollen, original REVEALS and optimized REVEALS data
		output: output/reconstructions/REVEALS_with_SD_forest_...
			output/reconstructions/pollen_forest_...

Technical specifications:
	- the code was written and run using R version 4.2.2 on Windows 10 x64
	- if problems occur while running the code please check your version