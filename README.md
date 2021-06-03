# brain_organoids

This project aims to identify enhancer-promoter important for lineage commitment. We therefore generated 5' single cell RNAseq at three timepoints of brain organoid development. First, we will use SCAFE to annotate the data. Second, we will use Harmony to integrate the datasets. In the next step we will use Seurat for clustering and cell type annotation, followed by pseudo-time analysis using Monocle to identify lineages. Using the DEG of the Monocle analysis, we correlated the expression values of enhancer and promoters along the pseudo-time. This correlation value represents how well the dynamics of the enhancer correlates with the promoter, with higher correlation values (ranging from 0 to 1) representing a more likely co-expression and therefore regulatory function of the enhancer to the pairing promoter. This initial pairing was further confirmed using Cicero, which is a tool predicting the co-accessibility of enhancers and promoters. 

Optional further analysis:
#^ use eQTL data 
#^ integrate HiC data for physical contact information
