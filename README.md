# ORBIT for *E. coli* : Kilobase-scale oligonucleotide recombineering at high throughput and high efficiency

Scott H Saunders, Ayesha M Ahmed

Green Center for Systems Biology and Lyda Hill Department of Bioinformatics

University of Texas Southwestern Medical Center, Dallas, Texas

## Abstract

Microbiology and synthetic biology depend on reverse genetic approaches to manipulate bacterial genomes; however, existing methods require molecular biology to generate genomic homology, suffer from low efficiency, and are not easily scaled to high throughput applications. To overcome these limitations, we developed a system for creating kilobase-scale genomic modifications that uses DNA oligonucleotides to direct the integration of a non-replicating plasmid. This method, Oligonucleotide Recombineering followed by Bxb-1 Integrase Targeting (ORBIT) was pioneered in *Mycobacteria*, and here we adapt and expand it for *E. coli*. Our redesigned plasmid toolkit achieved nearly 1000x higher efficiency than λ Red recombination and enabled precise, stable knockouts (\<134 kb) and integrations (\<11 kb) of various sizes. Additionally, we constructed multi-mutants (double and triple) in a single transformation, using orthogonal attachment sites. At high throughput, we used pools of targeting oligonucleotides to knock out nearly all known transcription factor and small RNA genes, yielding accurate, genome-wide, single mutant libraries. By counting genomic barcodes, we also show ORBIT libraries can scale to thousands of unique members (\>30k). This work demonstrates that ORBIT for *E. coli* is a flexible reverse genetic system that facilitates rapid construction of complex strains and readily scales to create sophisticated mutant libraries.

## Resources

-   [Targeting oligo design app](https://saunders-lab.shinyapps.io/ORBIT_TO_design_ecMG1655/)

-   [Saunders lab ORBIT help page](saunderslab.org/research/ec_orbit)

## Code repository

This github repository accompanies the manuscript and provides code notebooks that process and analyze raw data, turning it into preliminary figures in R. There are three main folders that are logically connected:

`/data -> /code -> /figures`

These code notebooks can be easily viewed as html documents. Links to these documents are listed below and are also available at the website version of this repository: <https://saunders-lab.github.io/ecoli_orbit/>

### `/data`

-   Low throughput experiments

-   High throughput experiments

-   seq_data

-   DNA maps

### `/code`

This folder contains the code notebooks used to process and analyze the data.

#### Main figures

-   Figure 1. ORBIT overview and proof of principle. (No computational content)

-   [Figure 2. ORBIT protocol optimization.](https://saunders-lab.github.io/ecoli_orbit/code/main_figs/fig_2_protocol_optimization.html)

-   [Figure 3. Benchmarking ORBIT.](https://saunders-lab.github.io/ecoli_orbit/code/main_figs/fig_3_gold_stds.html)

-   [Figure 4. Large deletions and insertions with ORBIT.](https://saunders-lab.github.io/ecoli_orbit/code/main_figs/fig_4_sizes.html)

-   [Figure 5. Orthogonal double mutants.](https://saunders-lab.github.io/ecoli_orbit/code/main_figs/fig_5_multi_orbit.html)

-   Figure 6. Markerless and scarless ORBIT strategies. (No computational content)

-   [Figure 7. Barcoded and multi-locus mutant libraries.](https://saunders-lab.github.io/ecoli_orbit/code/main_figs/fig_7_BC_oPool.html)

-   [Figure 8. High throughput knockout libraries for transcription factor and small RNA genes.](TBD) link TBD

#### Supplemental figures

-   Figure S1. Target oligo design principles and web app. (No computational content)

-   [Figure S2. Commercial oligo parameters and pricing.](https://saunders-lab.github.io/ecoli_orbit/code/supp_figs/supplement_protocol_optimization.html)

-   Figure S3. Troubleshooting attB only colonies. (No computational content)

-   Figure S4. Single gene deletion accuracy and stability. (No computational content)

-   [Figure S5. Helper plasmid efficiency and off target SNP mutation rate.](https://saunders-lab.github.io/ecoli_orbit/code/supp_figs/supp_mutation_rate.html)

-   Figure S6. Variable deletion size phenotypic accuracy. (No computational content)

-   Figure S7. Sanger sequence confirmation of *galK* deletions of various sizes. (No computational content)

-   Figure S8. Sanger sequence confirmation of *hisA, metA* and *leuD* deletions of various sizes. (No computational content)

-   Figure S9. Phenotypic accuracy of orthogonal and untargeted double mutants. (No computational content)

-   [Figure S10. Orthogonal att sites and triple mutant construction.]

-   [Figure S11. oPool library - Upstream and downstream reads for target loci.](https://saunders-lab.github.io/ecoli_orbit/code/supp_figs/mut_lib_grids.html)

-   Figure S12. Low abundance (Twist Biosciences) targeting oligo library processing. (No computational content)

-   [Figure S13. Short transcription factor library - Upstream and downstream reads for target loci.](https://saunders-lab.github.io/ecoli_orbit/code/supp_figs/mut_lib_grids.html)

-   [Figure S14. Long transcription factor deletion library - Upstream and downstream reads for target loci.](https://saunders-lab.github.io/ecoli_orbit/code/supp_figs/mut_lib_grids.html)

-   [Figure S15. Small RNA deletion library - Upstream and downstream reads for target loci.](https://saunders-lab.github.io/ecoli_orbit/code/supp_figs/mut_lib_grids.html)

-   [Figure S16. Explaining mutant library abundances.](https://saunders-lab.github.io/ecoli_orbit/code/supp_figs/supp_twist_lib.html)

-   [Figure S17. Correlations of oligo parameters and observed reads per target (perfect reads).](https://saunders-lab.github.io/ecoli_orbit/code/supp_figs/supp_twist_lib.html)

#### Sequencing data

##### galK BC

-   [galK random genomic barcodes - sequencing data analysis.](https://saunders-lab.github.io/ecoli_orbit/code/seq_data_processing/galK_BC_16N/galK_BC_16N_processing_analysis.html)

##### Twist oligo design

-   [Twist oligo pool design - Transcription factors](https://saunders-lab.github.io/ecoli_orbit/code/twist_oligo_design/twist_orbit_TF_deletion.html)

-   [Twist oligo pool design - Small RNAs](https://saunders-lab.github.io/ecoli_orbit/code/twist_oligo_design/twist_orbit_small_RNA_deletion.html)

-   [Twist oligo pool design - Primers and cloning scheme](https://saunders-lab.github.io/ecoli_orbit/code/twist_oligo_design/twist_orbit_cloning_scheme.html)

##### IDT oPool / Twist mutant library characterization

-   [Mutant libraries - Junction #1 sequencing data processing.](https://saunders-lab.github.io/ecoli_orbit/code/seq_data_processing/left_side_TO_libs_processing.html)

-   [Mutant libraries - Junction #2 sequencing data processing.](https://saunders-lab.github.io/ecoli_orbit/code/seq_data_processing/right_side_TO_libs_processing.html)

-   [Mutant libraries - Junction #1 analysis.](https://saunders-lab.github.io/ecoli_orbit/code/seq_data_processing/left_side_TO_libs_analysis.html)

-   [Mutant libraries - Junction #2 analysis.](https://saunders-lab.github.io/ecoli_orbit/code/seq_data_processing/right_side_TO_libs_analysis.html)

-   [Mutant libraries - Combined junction 1 & 2 analysis.](https://saunders-lab.github.io/ecoli_orbit/code/seq_data_processing/upstream_downstream_combined_analysis.html)

##### Twist TO abundance

-   [Twist TO abundance - Sequencing data processing.](https://saunders-lab.github.io/ecoli_orbit/code/seq_data_processing/TO_abundance/TO_abundance_processing.html)

-   [Twist TO abundance - Analysis.](https://saunders-lab.github.io/ecoli_orbit/code/seq_data_processing/TO_abundance/TO_abundance_analysis.html)

### `/figures`

This folder contains pdfs output from R containing preliminary figures. Final figures were made in illustrator.
