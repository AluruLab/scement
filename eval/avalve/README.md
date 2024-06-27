# Analysis of Human Aortic Valve dataset

Dataset is obtained from aortic valves of healthy and diseased individuals. 
Source paper : Kang Xu, Shangbo Xie,Yuming Huang,Tingwen Zhou, Ming Liu, Peng Zhu, Chunli Wang, Jiawei Shi, Fei Li,Frank W. Sellke and Nianguo Dong (2020) Cell-Type Transcriptome Atlas of Human Aortic Valves Reveal Cell Heterogeneity and Endothelial to Mesenchymal Transition Involved in Calcific Aortic Valve Disease.

The following table lists the scripts used for running the corresponding integration method:

| Method   | Script              |
| -------- | ------------------- |
| Combat   | avalve_combat.py    |
| FastMNN  | avalve_fastmnn.R    |
| Scanorma | avalve_scanorama.py |
| Seurat   | avalve_seurat.R     |
| SCEMENT  | avalve_scement.py   |

avalve_plots.py scripts generates plots for combat and scanorama.
avalve_data_pp.py is the script used for data pre-processing input data.
