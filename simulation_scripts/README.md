# Simulation script descriptions

This file describes the simulation scripts used to demonstrate the utility of PleistoDist, split up into simulation code, batch scripts for running simulations in parallel across multiple nodes, and plotting scripts. Each script has moderately detailed comments. If you heavily use these scripts to make your own code, consider citing the PleistoDist paper. Contact Ethan Gyllenhaal (egyllenhaal@unm.edu) if you have questions.

## Primary scripts

**msp_fiji_fluctuate.py:** Python script for dynamic msprime simulations, using python 3.7 and msprime v1. Takes in modified PleistoDist input and calculates parameters for IslandBioGeneflow and population size estimates between islands. Then simulates demography changing migration rates and population sizes each interval. Calculates F<sub>ST</sub> and nucleotide diversity between populations using simulated variants.

**msp_fiji_static.py:** Like above, but without dynamic change.

**samoa_dynamic.slim:** Eidos script for dynamic quantitative SLiM3 simulations. Uses a spatially explicit Wright-Fisher model with neutral mutations and individuals randomly spread across the Samoan archipelago. The map is generated from PleistoDist rasters converted to PNGs. Offspring can either undergo short- or long-distance dispersal. The latter calculates a random direction and distance, and succeeds if the individual never crosses water or crosses water and moves onto land (it stops at the first check where land is encountered). Individuals are assigned based on coordinates in Eidos, using X and Y coordinate checks instead of using a raster. After that, F<sub>ST</sub> and nucleotide diversity are calculated.

**samoa_static.slim:** Like above, but without dynamic change.

**solomons_colonization_nonWF.slim:** Similar to the dynamic Samoa script, but for the main Solomons archipelago (i.e., Bougainville through Makira), and using a non-WF model to allow for population size to change as the archipelago is colonized. Individuals start around one point on the west end of the archipelago. Note that this is made for qualitative visualization of the colonization process, and has no real output. It is meant to use SLiM3's excellent GUI to visualize the colonization history.

## Batch scripts

**run_samoa_dynamic.slurm:** Slurm batch script for using GNU parallel to run 100 replicates of each dispersal value for dynamic SLiM3 simulations, each on one core. Also initializes headers for output files. This can be scaled up to replicates * # dispersal values cores.

**run_samoa_static.slurm:** Like above, but for static simulations.

**run_fiji_msp.slurm:** Slurm batch script for using GNU parallel to run 100 replicates of each alpha and dispersal value for both types of msprime simulation, each on one core. Also initializes headers for output files. Note that this can't be scaled to core above the number of replicates without modifications.

## Plotting scripts

**plot_samoa_slim.R:** Script for making barplots for F<sub>ST</sub> and nucleotide diversity for select islands and island pairs. Pairs are chosen to demonstrate the effects of isolation and island size.

**plot_fiji_msp.R:** Script for making barplots for F<sub>ST</sub> and nucleotide diversity for the island pair.


## Spatial simulation GIFs

**glacial_colonization.gif:** GIF from SLiM GUI of colonzation starting during an interval near the last glacial maximum sea level.

![Glacial Colonization GIF](https://github.com/g33k5p34k/PleistoDistR/blob/main/simulation_scripts/glacial_colonization.gif)

**interglacial_colonization.gif:** GIF from SLiM GUI of colonzation starting during an interval near the last glacial maximum sea level.

![Interglacial Colonization GIF](https://github.com/g33k5p34k/PleistoDistR/blob/main/simulation_scripts/interglacial_colonization.gif)
