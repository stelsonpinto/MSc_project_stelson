MSc_project_stelson_Pinto


### MSc Project – Genetic links between autoimmune disease risk and infectious disease outcome
### aims: To identify shared susceptibility loci and candidate causal variants between systemic lupus erythematosus (SLE) and severe COVID-19 using Bayesian fine-mapping, Bayesian colocalisation and whole-blood cis-eQTL integration across European, East Asian and trans-ancestry datasets

### Project overview

This repository contains the code, metadata and analysis workflow used for my MSc Applied Bioinformatics dissertation at King's College London.

The project uses Bayesian fine-mapping (SuSiE), Bayesian colocalisation (coloc.susie), conditional analysis (GCTA-COJO) and whole-blood cis-eQTL integration to identify shared genetic susceptibility loci between systemic lupus erythematosus (SLE) and severe COVID-19 across European, East Asian and trans-ancestry datasets.

### Repository structure

- CODE/ – Analysis scripts.
- DATA/ – Metadata, reference files and sample lists.
- RESULTS/ – Analysis outputs and figures.
- ENVIRONMENT/ – Software environment and R session information.
- LOGS/ – Pipeline log files.
- docs/ – Workflow, software versions and data sources.

### Documentation

- `docs/Pipeline_workflow` – Analysis workflow.
- `docs/data_sources` – Summary of all datasets used.
- `docs/software_versions` – Software versions used.
- `docs/SAMPLE_sizes` – Sample sizes for all GWAS datasets.

### Data availability

Due to data licensing restrictions, the original GWAS, eQTL and reference genotype datasets are **not included** in this repository. Information on obtaining these datasets is provided in `docs/data_sources`.

### Reproducibility

All analyses were performed using the GRCh38 reference genome. Software versions are provided in `docs/software_versions`, and complete R package versions are available in `ENVIRONMENT/sessionInfo.txt`.
