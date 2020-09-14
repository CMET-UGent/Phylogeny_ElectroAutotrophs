# Phylogeny_ElectroAutotrophs

Underlying code and data for the data analysis of Prévoteau *et al.* (2020). 
Released under GPL-3. 

## File and folder structure

File/Folder | Contents 
------------|------------------
Source_data | Raw and processed datafiles used for downstream processing.
EPA.Rmd     | Code documention of the evolutionary placement algorithms and pre-processing of raw data. 
Processing.R| supporting code for pre-processing of data required for `EPA.Rmd`
[EPA.md](https://github.com/CMET-UGent/Phylogeny_ElectroAutotrophs/blob/master/EPA.md)      | Markdown rendering of `EPA.Rmd`
EPA_fileS/figure-html | Folder to store rendered figures for `EPA.md`
mothur_illumina_processing.log | log file containing commands and output for the mothur-based data-processing
Biocathode.files | mothur-compatible files file to indicate which files to merge (associated with genbank PRJNA641381)
bc.oligos | mothur-compatible oligos file with the primer-set used

## Dependencies

- [R](https://www.r-project.org) (>3.5.0.)
- R packages [ape](https://cran.r-project.org/web/packages/ape/index.html), [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html),  [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html), [DECIPHER](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html), [phyloseq](https://joey711.github.io/phyloseq/), [CMETNGS](https://github.com/CMET-UGent/CMETNGS),
[openxlsx](https://ycphs.github.io/openxlsx/index.html), [treeio](https://yulab-smu.github.io/treedata-book/) and the [tidyverse](https://www.tidyverse.org/)
- [Standard RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)
- [EPA-NG](https://github.com/Pbdas/epa-ng)
- [mothur](https://mothur.org/) (>1.42)
- BioEdit
- [genesis](https://groups.google.com/forum/#!searchin/raxml/visualize$20jplace$20itol%7Csort:date/raxml/IildAkMRI9Q/dBFEsUKlBQAJ)

## Publications cited and data obtained from/kindly provided by

- Eddie, Brian J., Zheng Wang, Anthony P. Malanoski, Richard J. Hall, Steve D. Oh, Cheryl Heiner, Baochuan Lin, and Sarah M. Strycharz-Glaven. "‘Candidatus Tenderia electrophaga', an uncultivated electroautotroph from a biocathode enrichment." International Journal of Systematic and Evolutionary Microbiology 66, no. 6 (2016): 2178-2185.
- Milner, Edward M., Dorin Popescu, Tom Curtis, Ian M. Head, Keith Scott, and H. Yu Eileen. "Microbial fuel cells with highly active aerobic biocathodes." Journal of Power Sources 324 (2016): 8-16.
- Xia, Xue, Yanmei Sun, Peng Liang, and Xia Huang. "Long-term effect of set potential on biocathodes in microbial fuel cells: electrochemical and phylogenetic characterization." Bioresource technology 120 (2012): 26-33.
- Strycharz-Glaven, Sarah M., Richard H. Glaven, Zheng Wang, Jing Zhou, Gary J. Vora, and Leonard M. Tender. "Electrochemical investigation of a microbial solar cell reveals a nonphotosynthetic biocathode catalyst." Applied and environmental microbiology 79, no. 13 (2013): 3933-3942.
- Rothballer, Michael, Matthieu Picot, Tina Sieper, Jan BA Arends, Michael Schmid, Anton Hartmann, Nico Boon, Cees JN Buisman, Frédéric Barrière, and David PBTB Strik. "Monophyletic group of unclassified γ-Proteobacteria dominates in mixed culture biofilm of high-performing oxygen reducing biocathode." Bioelectrochemistry 106 (2015): 167-176.
- Desmond-Le Quéméner, Elie, Mickaël Rimboud, Arnaud Bridier, Céline Madigou, Benjamin Erable, Alain Bergel, and Théodore Bouchez. "Biocathodes reducing oxygen at high potential select biofilms dominated by Ectothiorhodospiraceae populations harboring a specific association of genes." Bioresource Technology 214 (2016): 55-62.
- Liao, Chengmei, Jiali Wu, Lean Zhou, Tian Li, Qing Du, Jingkun An, Nan Li, and Xin Wang. "Optimal set of electrode potential enhances the toxicity response of biocathode to formaldehyde." Science of The Total Environment 644 (2018): 1485-1492.
- Prévoteau, Antonin, Peter Clauwaert, Frederiek-Maarten Kerckhof, and Korneel Rabaey. "Oxygen-reducing microbial cathodes monitoring toxic shocks in tap water." Biosensors and Bioelectronics 132 (2019): 115-121.

# Planned/not-yet executed deeper analyses of this data

Apart from the valid conclusions of Prevotéau *et al.* 2020, deeper data analysis
was planned but as-of-yet not executed. Interested community members are welcome
to do so, the authors will provide more info if needed.

- Full length guide tree with more representatives from Marinobacter-Chromatiaceae-Labrenzia (MCL) (e.g. from the SILVA LTP project)
- Building taxonomic trees for only V3-V4, V3 or V4 - this may not add much phylogenetic signal but no EPA is needed. The disadvantage is that a large number of reference publications don't have data in this region. 