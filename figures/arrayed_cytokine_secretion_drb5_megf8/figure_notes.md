
# Combined MEGF8 and DRB5 cytokine plot:

Objective - to compare DRB5 and MEGF8 cytokine data across 4 experiments, TCSL241, TCSL242, TCSL248, TCSL250. I want to compare them to CD3Z, CD28, and 4-1BB. 

## tcsl241 cytokines alone (NOT USING, just for reference)
	https://benchling.com/s/etr-XVHoz7Y7e6Qnwdggl9sw?m=slm-R2D6cVHLrcyC1eZr9CT0

	- individual amounts, in python 
		/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/240202_tcsl241_cytokines
		raw data is: 
		/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/240202_tcsl241_cytokines/240202_TCSL241_Cytokines_Default Unmixed Worksheet.CSV
		previous python plots, which describe the data format: 
		/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/240202_tcsl241_cytokines/cytokine_table.py


## tcsl242 cytokines alone (NOT USING, just for reference)

	- final cytokines at end of 5 Repstims:
	/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/240328_tcsl242_cytokine/plot_cytokines.Rmd

## combined tcsl241 and tcsl242 data:

	raw csv: /Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/240311_tcsl_242_cytokine/240311_cytokine_both_donors.csv
	R file: /Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/240311_tcsl_242_cytokine/cytokines.Rmd


## tcsl247:
	extra cells from tcsl245 in vivo, as a 3rd replicate, contains DRB5 and DRB528C, which is called DRB5.CD28 in this experiment. does not contain MEGF8. 
	This experiment has two donors as well as a separate replicate of donor 1 called 1V which are the left over sorted in vivo cells. See the metadata file. I only have zeta for the extra in vivo cell set. 

	- raw data from gating on the aurora
	/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20240528_tcsl247_invitro_cytokine/TCSL247 Cytokine Unmixed.csv
	- plate map
	/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20240528_tcsl247_invitro_cytokine/plate_map.csv
	- car name map
	/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20240528_tcsl247_invitro_cytokine/sample_metadata.csv


## tcsl248: 

	metadata and many files in:
	/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241112_tcsl248_cytokine/

		- raw cytokine data
		/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241112_tcsl248_cytokine/241112_tcsl248_cytokine_TCSL248 Cytokine Unmixed.CSV

		- flowjo gated data 
		/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241112_tcsl248_cytokine/populations.csv

		- Rmarkdown plot for this experiment only
		/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241112_tcsl248_cytokine/cytokine_analysis.Rmd

		- plate map for well data
		/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241112_tcsl248_cytokine/plate_map.csv

		- CAR map for labelling samples
		/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241112_tcsl248_cytokine/car_map.csv

## tcsl250: 

	metadata and many files in:
	/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241213.tcsl250 cytokines/

		- raw cytokine data 
		/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241213.tcsl250 cytokines/2024.12.13 TCSL250 Cytokines_TCSL247 Cytokine Unmixed.CSV

		- flowjo gated data 
		/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241213.tcsl250 cytokines/populations.csv

		- Rmarkdown plot for this experiment only
		/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241213.tcsl250 cytokines/cytokine_analysis.Rmd

		- plate map for well data
		/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241213.tcsl250 cytokines/plate_map.csv

		- CAR map for labelling samples
		/Users/dbg/Library/CloudStorage/Box-Box/tcsl/flow/20241213.tcsl250 cytokines/car_map.csv


Each of these has a slightly different format, with 248 and 250 being more similar. We need to collect and merge all the data from these experiments and make an integrated table and plot of the cytokine secretion across all 4 donors for CD3Z, CD28, MEGF8, and DRB5. 
