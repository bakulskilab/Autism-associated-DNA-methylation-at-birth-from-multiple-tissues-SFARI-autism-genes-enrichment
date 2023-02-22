# Autism-Associated DNA Methylation at Birth From Multiple Tissues Is Enriched for Autism Genes in the Early Autism Risk Longitudinal Investigation

### Bakulski KM, Dou JF, Feinberg JI, Aung MT, Ladd-Acosta C, Volk HE, Newschaffer CJ, Croen LA, Hertz-Picciotto I, Levy SE, Landa R, Feinberg AP, and Fallin MD

## Citation Information

Bakulski, K. M., Dou, J. F., Feinberg, J. I., Aung, M. T., Ladd-Acosta, C., Volk, H. E., Newschaffer, C. J., Croen, L. A., Hertz-Picciotto, I., Levy, S. E., Landa, R., Feinberg, A. P., & Fallin, M. D. (2021). Autism-Associated DNA Methylation at Birth From Multiple Tissues Is Enriched for Autism Genes in the Early Autism Risk Longitudinal Investigation. Frontiers in molecular neuroscience, 14, 775390. https://doi.org/10.3389/fnmol.2021.775390. PMID: 34899183

This Github repository contains the data management and analytic scripts to produce the following manuscript: [Autism-Associated DNA Methylation at Birth From Multiple Tissues Is Enriched for Autism Genes in the Early Autism Risk Longitudinal Investigation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8655859/)


## Abstract

Background: Pregnancy measures of DNA methylation, an epigenetic mark, may be associated with autism spectrum disorder (ASD) development in children. Few ASD studies have considered prospective designs with DNA methylation measured in multiple tissues and tested overlap with ASD genetic risk loci.

Objectives: To estimate associations between DNA methylation in maternal blood, cord blood, and placenta and later diagnosis of ASD, and to evaluate enrichment of ASD-associated DNA methylation for known ASD-associated genes.

Methods: In the Early Autism Risk Longitudinal Investigation (EARLI), an ASD-enriched risk birth cohort, genome-scale maternal blood (early n = 140 and late n = 75 pregnancy), infant cord blood (n = 133), and placenta (maternal n = 106 and fetal n = 107 compartments) DNA methylation was assessed on the Illumina 450k HumanMethylation array and compared to ASD diagnosis at 36 months of age. Differences in site-specific and global methylation were tested with ASD, as well as enrichment of single site associations for ASD risk genes (n = 881) from the Simons Foundation Autism Research Initiative (SFARI) database.

Results: No individual DNA methylation site was associated with ASD at genome-wide significance, however, individual DNA methylation sites nominally associated with ASD (P < 0.05) in each tissue were highly enriched for SFARI genes (cord blood P = 7.9 × 10<sup>–29</sup>, maternal blood early pregnancy P = 6.1 × 10<sup>–27</sup>, maternal blood late pregnancy P = 2.8 × 10<sup>–16</sup>, maternal placenta P = 5.6 × 10<sup>–15</sup>, fetal placenta P = 1.3 × 10<sup>–20</sup>). DNA methylation sites nominally associated with ASD across all five tissues overlapped at 144 (29.5%) SFARI genes.

Conclusion: DNA methylation sites nominally associated with later ASD diagnosis in multiple tissues were enriched for ASD risk genes. Our multi-tissue study demonstrates the utility of examining DNA methylation prior to ASD diagnosis.

## Funding
Funding for the EARLI study was provided by NIH (R01ES016443, PI: CN) and Autism Speaks. The DNA methylation measures were funded by NIH (R01ES017646, PIs: MF and AF). MF, JD, and KB were supported by grants from the National Institute of Environmental Health Sciences (R01 ES025531; R01 ES025574; and P30 ES017885). KB was supported by a grant from the National Institute of Aging (R01 AG055406).

## Script Files
ASD_sens_functions.R: filtering CpG sites and preparing for sensiticity analysis

ASDsensitivity.R: code for sensitivity analysis

