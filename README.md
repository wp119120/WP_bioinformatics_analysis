# Integrated Transcriptomic, Machine Learning, and Molecular Dynamics Analysis Pipeline for PFOA-Induced Pulmonary Fibrosis

## Overview

This repository contains the complete analysis pipeline, analytical scripts, and molecular simulation configurations used in our study investigating the pathogenic mechanisms, diagnostic biomarkers, and structural binding properties of Perfluorooctanoic Acid (PFOA) in pulmonary fibrosis. The computational framework seamlessly integrates multicenter transcriptomics (GEO batch correction), 127 ensemble machine learning algorithms, immune infiltration profiling, molecular docking, and 100 ns all-atom molecular dynamics (MD) simulations.

---

## 🧬 Study Background

Pulmonary fibrosis is a progressive, irreversible respiratory disorder characterized by alveolar structural destruction and excessive extracellular matrix deposition. Perfluorooctanoic acid (PFOA), a pervasive and non-biodegradable environmental toxicant, has been increasingly recognized as a hazardous contributor to respiratory inflammation and tissue remodeling. This repository provides an end-to-end computational biology workflow identifying pivotal diagnostic biomarkers (LCN2, TYMS, BCAT2, NQO1, PRDX6, CLPP) and elucidating spontaneous physical-chemical target interactions via structural biology.

---

## 🔬 Study Design and Computational Workflow

```text
               [ Multicenter Microarray Cohorts (GEO) ]
                                   │
                                   ▼
          [ ComBat Batch Correction & limma DEG Identification ]
                                   │
                                   ▼
         [ PFOA Target Prediction & Adverse Target Intersection ]
                                   │
                                   ▼
     [ 127 Ensemble Machine Learning Models & SHAP Interpretability ]
                                   │
                                   ▼
         [ Diagnostic Nomogram & CIBERSORT Immune Landscape ]
                                   │
                                   ▼
     [ AutoDock Vina Docking & 100 ns GROMACS MD Simulations ]
## 📁 Repository Structure

```text
├── 01_Bioinformatics_and_Machine_Learning/
│   ├── 01_Batch_Correction_and_PCA.R          # SVA ComBat harmonization & PCA
│   ├── 02_Differential_Expression.R           # limma DEG screening, Volcano & Heatmap
│   ├── 03_Target_Venn_Intersection.R          # PFOA targets vs DEGs intersection
│   ├── 04_GO_Enrichment_Analysis.R            # Functional annotation & Circos visualization
│   ├── 05_ML_Data_Preprocessing.R             # Train/Test cohort partitioning
│   ├── 06_Run_127_Machine_Learning.R          # 127 algorithmic pipeline benchmarking
│   ├── 07_SHAP_Interpretability.R             # SHAP beeswarm & feature dependence
│   ├── 08_Core_Genes_ROC_and_Boxplot.R        # Diagnostic ROC curves & expression boxplots
│   ├── 09_Nomogram_and_DCA_Analysis.R         # Nomogram calibration & Decision Curve Analysis
│   ├── 10_CIBERSORT_Immune_Infiltration.R     # SVR immune deconvolution
│   ├── 11_Immune_Infiltration_Visualization.R # Relative fraction stacked barplots
│   ├── 12_Immune_Differential_Analysis.R      # Wilcoxon differential tests across 22 phenotypes
│   └── 13_Immune_Correlation_and_Network.R    # Biomarker-immune correlation & LinkET network
├── 02_Molecular_Docking/
│   ├── 01PFOA-TYMS/                           # Receptor PDBQT, ligand, grid parameters & docked poses
│   ├── 02PFOA-LCN2/                           # Receptor PDBQT, ligand, grid parameters & docked poses
│   ├── 03PFOA-NQO1/                           # Receptor PDBQT, ligand, grid parameters & docked poses
│   ├── 04PFOA-PRDX6/                          # Receptor PDBQT, ligand, grid parameters & docked poses
│   ├── 05PFOA-CLPP/                           # Receptor PDBQT, ligand, grid parameters & docked poses
│   └── 06PFOA-BCAT2/                          # Receptor PDBQT, ligand, grid parameters & docked poses
├── 03_Molecular_Dynamics_Simulation/
│   ├── 01PFOA-TYMS/                           # Coordinates (GRO), parameters (MDP), & trajectory (XTC)
│   ├── 02PFOA-LCN2/                           # Coordinates (GRO), parameters (MDP), & trajectory (XTC)
│   ├── 03PFOA-NQO1/                           # Coordinates (GRO), parameters (MDP), & trajectory (XTC)
│   ├── 04PFOA-PRDX6/                          # Coordinates (GRO), parameters (MDP), & trajectory (XTC)
│   ├── 05PFOA-CLPP/                           # Coordinates (GRO), parameters (MDP), & trajectory (XTC)
│   └── 06PFOA-BCAT2/                          # Coordinates (GRO), parameters (MDP), & trajectory (XTC)
└── README.md
```