# SLiM simulations to test the mutation area idea that Moi outlined


The SLiM simulation configuration file in ```Slim_configs/simulationDevelopment_extinctionColonisationTemplate.slim``` contains a model of a 2D stepping stone meta-population. Local adaptation to environmental heterogeneity is incorporated under an arbitrart genetic architechture of adaptation.

At a particular point in the simulation, we model habitat destruction by making it impossible for individuals to survive in certain locations. The idea is to test how habitat loss will impact genetic variation and, particularly, adaptability. At this stage, additive genetic variance is used to model adaptability.

The Simulation is run as a typical SLiM simulation. The simulation outputs a VCF file for each population at specified intervals. These are then processed and compressed to reduce file size.

Here is the pipeline I used to:
1. Run the simulation
2. Process the output
3. Clean up

```sh
mkdir test_run
cd test_run
slim ../Slim_configs/simulationDevelopment_extinctionColonisationTemplate.slim
mkdir populationVCFs
mv *vcf populationVCFs/

parallel "python ../bin/summariseSLiMulation.py populationVCFs/{} {}.summary.csv" ::: $(ls populationVCFs/ | cut -f1 -d'.' |uniq)

mkdir alleleFreq_summaries  phenotype_summaries

mv alleleFreqs* alleleFreq_summaries/

mv phenotypes.* phenotype_summaries/

gzip populationVCFs/*

cd ../
```

Now we have a set of files that contain information on the phenotypic, additive genetic and genetic diversity over time.
