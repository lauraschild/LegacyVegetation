Code to reconstruct vegetation from a subset of pollen records following the LegacyVegetation
dataset. In this repository parts of the analysis can be replicated using a subset of three records.

The script main_complete.R includes all of the analyses and takes approximately 10 minutes to run all
scripts in the necessary order. It is recommended to run this first as it also takes care of packages
that need to be installed.

The following scripts are included in this analysis:
	- main_REVEALS.R
		reconstructs vegetation cover from three pollen records using the REVEALS model and the
		original RPP estimates
		output: output/PANGAEA/Europe_original_RPP.csv
	- main_optimization.R
		optimizes RPP for one taxon (ten in manuscript) to fit modern reconstructed forest cover to 
		remote sensing forest cover
		output: output/optimized_RPP_Europe.csv
			output/figures/valid_Europe.png
	- main_opti_REVEALS.R
		uses the optimized RPP to run REVEALS on original pollen record
		output: output/PANGAEA/optimized_REVEALS_Europe.csv
	- main_reconstruct_forest.R
		reconstruct forests cover from Pollen, original REVEALS and optimized REVEALS data
		output: output/PANGAEA/composition_forest...
			output/PANAEA/forest...

The results of the optimization will not be comparable to the results shown in the manuscript, as only 
three records are being used to optimize just one taxon.

Technical specifications:
	- the code was written and run using R version 4.2.2 on Windows 10 x64
	- if problems occur while running the code please check your version