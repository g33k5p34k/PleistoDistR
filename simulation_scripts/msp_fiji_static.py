'''
Made by Ethan Gyllenhaal (egyllenhaal@unm.edu)
Last updated 12 April 2022

Script for running a single msprime simulation for static divergence between two islands (designed for parallel)
Takes in parameters for infile path, output path, alpha (# propagules per unit^2), and mean dispersal distance
Infile is comma delimited with header Interval,Distance,Kadavu_target,Viti_target,Kadavu_size,Viti_size
Outputs nucleotide diversity for each island and FST

calcParameters calculates gene flow and population size per interval given input file.
run_msprime is the method for running the simulation and calculating summary statistics
main is simply the driver for everything
'''

import os, sys, math, msprime as msp, numpy as np, re, allel, argparse, pandas as pd

# global for how many generations there are per interval
# assuming a generation time of 2 years and 5k years/interval
generation_interval = 2500
# global for population density per unit^2 (meter here)
density = 0.00001

def main():
    
    # set up parser and arguments with ArgParse
    parse = argparse.ArgumentParser(description = "Get simulation parameters")
    parse.add_argument("-i", "--infile", type=str, help="Path to island information file")
    parse.add_argument("-o", "--output", type=str, help="Path to output file")
    parse.add_argument("-a", "--alpha", type=float, help="Value of alpha, propagules/unit^2")
    parse.add_argument("-d", "--dispersal", type=float, help="Value of mean dispersal distance in default units")
    args = parse.parse_args()
   
    # assign arguments to variables
    infile, outfile, alpha, dispersal = args.infile, args.output, args.alpha, args.dispersal
    
    # calculate the population sizes and migration rates over time
    # uses calcParameters method below, used as input for msprime
    processed_data = calcParameters(infile, alpha, dispersal, density)

    # amount of generations to run the simulation for
    time = 200000

    # calls msprime, assigns results to variables, and writes those to output
    fst, div1, div2 = run_msprime(processed_data, time)
    output = open(outfile, "a")
    output.write(str(fst)+"\t"+str(div1)+"\t"+str(div2)+"\n")
    output.close()
    
def run_msprime(params, time):
    # set number of samples
    samples=50
    # determine the interval the split occurs at
    split_int = time//generation_interval

    # pairing indexes and rows
    params = params.reset_index()

    # Set up populations, with pop1 as Viti Levu and pop2 as Kadavu
    demography = msp.Demography()
    demography.add_population(name="pop1", initial_size=params["Viti_pop"][0])
    demography.add_population(name="pop2", initial_size=params["Kad_pop"][0])
    demography.add_population(name="anc_pop12", initial_size = params["Viti_pop"][split_int]+params["Kad_pop"][split_int])
    demography.add_population_split(time=time, derived=["pop1","pop2"], ancestral="anc_pop12")
    
    # no dynamism! so just set migration rates straight up
    demography.set_migration_rate(source="pop2", dest="pop1", rate=params["K2V"][0]/params["Kad_pop"][0])
    demography.set_migration_rate(source="pop1", dest="pop2", rate=params["V2K"][0]/params["Viti_pop"][0])

    # sort the demographic events
    demography.sort_events()

    # Run the simulation to get tree sequences
    trees=msp.sim_ancestry(samples={"pop1":samples, "pop2":samples}, 
                           demography=demography, 
                           recombination_rate=1e-8,
                           sequence_length=1e7)
    
    # Add mutations to treeseqs
    mutation=msp.sim_mutations(trees, rate=2.3e-9)
    
    # get haplotypes and genotypes from simulation 
    haplotypes = np.array(mutation.genotype_matrix())
    positions = np.array([s.position for s in trees.sites()])
    genotypes = allel.HaplotypeArray(haplotypes).to_genotypes(ploidy=2)
    
    # calculate fst
    fst = allel.stats.fst.average_weir_cockerham_fst(genotypes,[list(range(0,int(samples))),list(range(int(samples),samples*2))],10)[0]
    
    # calculate nucleotide diversity
    geno1 = genotypes.take(range(0,int(samples)),axis=1)
    geno2 = genotypes.take(range(int(samples),int(samples*2)), axis=1)
    acount1 = geno1.count_alleles()
    acount2 = geno2.count_alleles()
    div1 = allel.sequence_diversity(range(1,len(acount1)),acount1)
    div2 = allel.sequence_diversity(range(1,len(acount2)),acount2)
    
    # return the summary statistics
    return(fst, div1, div2)

def calcParameters(infile, alpha, mean_disp, dens):
    # read in PleistoDist output as a dataframe
    dframe = pd.read_csv(infile, index_col=0)
    # Calculate population densities
    dframe["Viti_pop"] = dframe["Viti_size"] * density
    dframe["Kad_pop"] = dframe["Kadavu_size"] * density
    # BACKWARD TIME Viti -> Kadavu, i.e. forward Kadavu -> Viti
    dframe["V2K"] = (alpha * dframe["Kadavu_size"]) * ((dframe["Viti_target"] * np.exp(-1 * (1/mean_disp) * dframe["Distance"]))/
                                                       (2 * math.pi * dframe["Distance"]))
    # BACKWARD TIME Kadavu -> Viti, i.e. forward Viti -> Kadavu
    dframe["K2V"] = (alpha * dframe["Viti_size"]) * ((dframe["Kadavu_target"] * np.exp(-1 * (1/mean_disp) * dframe["Distance"]))/
                                                     (2 * math.pi * dframe["Distance"]))

    return dframe

if __name__ == '__main__':
    main()
