# Analysis of A. thaliana root dataset

Dataset is obtained from A. thaliana roots under various conditions. 
Source publications :
1.  E-GEOD-152766: Shahan R, Hsu C, Nolan TM, Cole BJ, Taylor IW et al. (2020) A single cell Arabidopsisroot atlas reveals developmental trajectories in wild type and cell identity mutants.
2.  E-GEOD-121619: Jean-Baptiste K, McFaline-Figueroa JL, Alexandre CM, Dorrity MW, Saunders L et al. (2019) Dynamics of Gene Expression in Single Root Cells of Arabidopsis thaliana.
3.  E-GEOD-123013: Ryu KH, Huang L, Kang HM, Schiefelbein J. (2019) Single-Cell RNA Sequencing Resolves Molecular Relationships Among Individual Plant Cells.
4.  E-GEOD-158761: Gala HP, Lanctot A, Jean-Baptiste K, Guiziou S, Chu JC et al. (2020) A single cell view of the transcriptome during lateral root initiation in Arabidopsis thaliana.

The following table lists the scripts used for running the corresponding integration method:

| Method            | Script              |
| ----------------- | ------------------- |
| Combat & SCEMENT  | ath_combat.py    |
| FastMNN           | ath_fastmnn.R    |
| Scanorma          | ath_scanorama.py |
| Seurat             | ath_seurat.R     |

ath_data_pp.py is the script used for data pre-processing input data.
