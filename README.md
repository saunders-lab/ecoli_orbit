# ORBIT for *E. coli* : Kilobase-scale oligonucleotide recombineering at high throughput and high efficiency

Scott H Saunders, Ayesha M Ahmed

Green Center for Systems Biology and Lyda Hill Department of Bioinformatics

University of Texas Southwestern Medical Center, Dallas, Texas

## Abstract

Microbiology and synthetic biology depend on reverse genetic approaches to manipulate bacterial genomes; however, existing methods require molecular biology to generate genomic homology, suffer from low efficiency, and are not easily scaled to high throughput applications. To overcome these limitations, we developed a system for creating kilobase-scale genomic modifications that uses DNA oligonucleotides to direct the integration of a non-replicating plasmid. This method, Oligonucleotide Recombineering followed by Bxb-1 Integrase Targeting (ORBIT) was pioneered in *Mycobacteria*, and here we adapt and expand it for *E. coli*. Our redesigned plasmid toolkit achieved nearly 1000x higher efficiency than λ Red recombination and enabled precise, stable knockouts (\<134 kb) and integrations (\<11 kb) of various sizes. Additionally, we constructed multi-mutants (double and triple) in a single transformation, using orthogonal attachment sites. At high throughput, we used pools of targeting oligonucleotides to knock out nearly all known transcription factor and small RNA genes, yielding accurate, genome-wide, single mutant libraries. By counting genomic barcodes, we also show ORBIT libraries can scale to thousands of unique members (\>30k). This work demonstrates that ORBIT for *E. coli* is a flexible reverse genetic system that facilitates rapid construction of complex strains and readily scales to create sophisticated mutant libraries.

## Resources

-   [Targeting oligo design app](https://saunders-lab.shinyapps.io/ORBIT_TO_design_ecMG1655/)

-   [Saunders lab ORBIT help page](https://saunders-lab.github.io/ecoli_orbit/)

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

-   [Figure 8. High throughput knockout libraries for transcription factor and small RNA genes.](TBD)

#### Supplemental figures

-   Figure S1.

### `/figures`

This folder contains pdfs output from R containing preliminary figures. Final figures were made in illustrator.
