import vcf
import sys, glob
import numpy as np
import pandas as pd
from functools import reduce

# Read the list of optima from the BC map
optima = [int(i.strip()) for i in open("/home/booker/work/MutationAreaSlimulations/Slim_configs/BC_Map_environments.14x14.txt","r").readlines()]

# initiate a dict to store phenotypic effects
s_by_pos_dict = {}

# initiate lists to store allele frequency and phenotypic information
allele_freq_dataFrames = []
phenotype_dataFrames = []

# Loop over all the VCF files that the simulation generated
for vcf_file in glob.glob(sys.argv[1] + "*vcf"):

# Get the population id from the file name
    this_pop = int(vcf_file.split("/")[-1].split(".")[2])

# Get the population optima from the list
    this_optimum = optima[this_pop-1]

# initiate a VCF Reader
    vcf_reader = vcf.Reader( open( vcf_file, "r" ) )

# initiate containters for phenotypes and allele frequencies
    phen_list = []
    alleleFreq_mat = []

# Iterate over each record in the VCF files
    for record in vcf_reader:
# grab the phenotpyic effect (this is labelled "s" by default in SLiM)
        phenotypicEffect = record.INFO["S"]
# add the phenotypic effect to the relevant dict
        s_by_pos_dict[record.POS] = phenotypicEffect[0]
# Widen the list of genotypes to get a list of haplotypes
        genotypes_wide = [ sample["GT"].split("|") for sample in record.samples]
# Re-flatten the list of lists
        genotypes = np.array( [int(hap) for sample in genotypes_wide for hap in sample] )

# make a temporary object that corresponds to the allelic contribution to the phenotype for each haplotype
        temp_phenotypes = genotypes*phenotypicEffect[0]

# add the allele freuqencies to the container dict
        alleleFreq_mat.append( [record.POS, genotypes.mean() ] )

# If a variant is neutral, don't bother calculating phenotypic contribution - it has none
        if phenotypicEffect[0] == 0.:
            #Neutral
            pass

        else:
# If we're dealing with a phenotype affecting variant, calculate the contribution to the individuals phenotype from this locus
            phen_list.append(  [ sum(temp_phenotypes[i:i+2]) for i in range(0, len(temp_phenotypes), 2)] )

# Make a pandas DataFrame from the phenotypes from this population
    phen_DF = pd.DataFrame(phen_list, columns = [str(n) for n in range(int( len(temp_phenotypes)/2))]).sum(axis = 0)

# Calculate the mean and variance of the phenotypes, add the population optimum and id to the dataframe for reference later
    phen_DF["meanPhen"] = phen_DF.mean(axis = 0)
    phen_DF["phenVar"] = phen_DF.var(axis = 0)
    phen_DF["optimum"] = this_optimum
    phen_DF["pop"] = this_pop

# Make a temporary DF from the allele frequnecy matrix for this population
    temp_DF = pd.DataFrame(alleleFreq_mat, columns = ["POS","p_"+str(this_pop)])
# Add the temporary DFs to the respective containers
    allele_freq_dataFrames.append( temp_DF )
    phenotype_dataFrames.append( pd.DataFrame(phen_DF).transpose())

# Make a single DF from phenotype information from each Pop
phenotype_df_merged = pd.concat(phenotype_dataFrames)
phenotype_df_merged = phenotype_df_merged.reindex( sorted(phenotype_df_merged.columns), axis = 1)

# Save to file
phenotype_df_merged.to_csv("phenotypes."+sys.argv[2], index = False)

# Merge all the allele freuqnecy DFs (sorry for the horrible Pandas one-liner)
allele_df_merged = reduce( lambda left, right: pd.merge(left, right, on = ["POS"], how = "outer"), allele_freq_dataFrames).fillna(0.0)

# Add the phenotypic effects to the Dataframe
allele_df_merged["s"] = allele_df_merged["POS"].map(s_by_pos_dict)

# Save to File
allele_df_merged.to_csv("alleleFreqs."+sys.argv[2], index = False)
