# Set of scripts associated with processing the data from the ForestGEO plots

*Note that all data is stored on the google drive itself, not on the GitHub

Scripts included (Run in outputed order)
00_Functions.R - file of commonly used functions throughout this analysis,
	sourced directly in other scripts as needed

AA_LocationManager.R - Script caleld by all other scripts to ensure consistency in location across all the scripts
	for both the data as well as other scripts.
 
00_PackageInstall.R - Script to ensure that all packags are present

00_DataInitialClean_XXXX.R - script dedicated for each site (e.g. _HVDF is for Harvard forest)
	to cleaning up and manipulating the data in a format that is specific to it. Many of
	the sites have different quirks for them. Overall layout is similar across all. Site names
	are as follows:
		HVDF - Harvard Forest
		SERC - Smithsonian Environmental Research center
		SCBI - Smithsonian Conservation Biology Institute
		UMBC - University of Maryland Baltimore County (data NOT used due to lack of forested area)
		ZOFN - Zofin (data not used due to restriction to US plots)

00_RemformatStdz.R - Restandardizes the data (done after initial run through due to discrepencies of dealign with the data)

01_SpatialFilter.R Goes through the sites and reworks data due to large portions of the plot that were not measured
	due to Ruger's algorithm requiring measurements of all stems and partitioning parts of the plots into
	subplots, any gaps are problematic and must be accounted for. This will remove some of the stems for the
	entire analysis.

01A_TALLOAllometries.R - generate allometries from TALLO data to get estimates of height from DBH

01B_CanopyAssignments.R - using estimated heights, assign trees to individual canopy layers

02_GetValues.R - retrieving specific growth, mortality and size information for each species at the site.

03_GenrtSTANData.R - reformat the data into a list for stan to actually use, outputs

04_STAN_AnalysisRunsOnly.R - actually runs the STAN models, does some bare bones model checks (see notes in script) during the process, outputs the 
	posterior distribution of the results.

05_AddFunctionalTraitData.R - Reads in the Functional trait data, gets it ready so it can be ready for subsequent scripts. Assumes that functional
	trait data is at the species leve rather than site or individual specific.

05_wPCA.R - running of the weighted PCA, generates plots and results along with some diagnostics. Can be configured to output results with or
	without functional trait data 


