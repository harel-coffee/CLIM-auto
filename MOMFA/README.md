# CLIM/MOMFA
Multiobjective Metabolic Flux Analysis Algorithm

genesymbol2HGNC.m: MATLAB function to convert gene symbols to HGNC gene IDs (identifier for genes-to-reaction relations in Recon 2.2).


# Input
modelSpecs_NCI60: Boundary or extracellular fluxes estimated from CORE dataset for ovarian cancer cell lines (NCI 60 panel). Input is formatted as input for COBRA toolbox for Recon 2.2 model

Metabolic models: Context-specific metabolic models reconstructed using gene expression data from TCGA are available in /CLIM/ModelRecon/OutputModels

Gene_Names.txt: Gene symbols and corresponding HGNC IDs
