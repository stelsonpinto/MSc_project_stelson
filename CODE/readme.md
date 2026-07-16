### CODE

This directory contains the scripts used to perform the computational analyses described in the dissertation.

Scripts are executed in the order.

Script | Description 

01_extract_loci.sh | Extract lead SNPs and genomic regions for downstream analyses.

02_extract_sle_loci.R | Extract SLE GWAS variants within each locus. 

02_extract_covid_loci.R | Extract severe COVID-19 GWAS variants within each locus. 

03_extracting_eQTLGen_region.sh | Extract corresponding cis-eQTLGen regions. 

04_liftover_eqtl.R | Convert eQTL coordinates to GRCh38 where required. 

05_align_GWAS_eQTL.R | Harmonise GWAS and eQTL datasets and align alleles. 

06_create_ld_reference.sh | Create ancestry-specific LD reference datasets. 

07_compute_LD.sh | Generate LD matrices for each analysed locus. 

FINE_MAPPING_LOOP.R | Run SuSiE fine-mapping across all loci. 

08_summarise_coloc_results.R | Summarise Bayesian colocalisation results. 

09_coloc_susie.R | Perform SuSiE fine-mapping and Bayesian colocalisation. 

METAL FOR COVID META | Using METAL software to merge covid GWAS summary statistics of european and eas asian ancestry. 


Detailed execution order and analysis workflow are provided in `docs/Pipeline_workflow`.
