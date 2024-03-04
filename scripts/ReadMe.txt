Overview of Scripts:

List of Scripts:
00_Functions
- Common functions used throughout the scripts, especially in the first portion
that involves a lot of clean-up

00_PackageInstall
- List of packages used throughout the entire project. Note that this does NOT 
include any packages for stan/bayesian plotting, since STAN can be rather finicky
about some other stuff unfortunately.

00_DataInitialClean_XXXX (Where XXXX is a four letter code)
- Code specifically designed for each site to clen up the data into a proper format.
Individual sites have some small differences in terms of data format, so this works to
standardize everyinthg.

00_ReformatStdz
- Code to reformat the data to a new format that allows for 

01A_TALLOAllometries
- Takes the data from the TALLO publication of Jucker et al. 2022,
 Global Change Biology, doi: 10.1111/gcb.16302 (doi to specific dataset are in the
listed in the 2022 publication)