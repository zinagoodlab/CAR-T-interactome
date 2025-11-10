# CAR-T-interactome

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-pending-blue.svg)](https://doi.org/pending)

Single-cell interactome analysis reveals CCL8⁺CCL13⁺ tumor-associated macrophages as drivers of CAR T cell resistance in large B cell lymphoma. This repository contains computational pipelines for identifying resistance mechanisms through on-treatment single-cell profiling and cross-cohort validation.
<br>
<br>
<img width="5244" height="5243" alt="Cancer Cell Graphical Abstract" src="https://github.com/user-attachments/assets/b1397334-068d-4ac9-9178-0c24b7fc129c" />
<br>
<br>
We performed the first on-treatment single-cell characterization of the tumor microenvironment following CAR T cell therapy in large B cell lymphoma (LBCL). Through comprehensive single-cell multi-omic profiling and interactome analysis, we identified CCL8⁺CCL13⁺ tumor-associated macrophages (TAMs) as key mediators of CAR T resistance.
<br>
<br>
<img width="3106" height="1125" alt="image" src="https://github.com/user-attachments/assets/8e1e1470-75ba-43e1-ba87-81d72581fb24" />
<br>
<br>
**Key Findings**
- CCL8⁺CCL13⁺ TAMs are highly proliferative and immunosuppressive, predicting poor outcomes across multiple treatment contexts
- IRF4/IRF8-driven dysfunction in CAR T cells operates independently of canonical TOX-mediated exhaustion
- CCL18⁺ TAMs promote immunosuppression via CCR8⁺ regulatory T cells
- Resistance mechanisms are conserved between CAR T therapy and chemoimmunotherapy
  
**Dataset**
- Discovery cohort: 18 LBCL patients + 1 B-ALL patient receiving CAR T therapy
_ Validation cohorts: 785 baseline samples from 6 independent LBCL cohorts
_ Single-cell data: 31,808 high-quality cells from on-treatment tumor biopsies
- Technologies: scRNA-seq, scTCR-seq, CITE-seq using 10X Genomics platform
  
**Interactome Analysis Pipeline**
- REMI
- CellPhoneDB
- CellChat
  
**Validation**
- CIBERSORTx deconvolution of bulk RNA-seq datasets
- Cross-cohort validation across CAR T and chemoimmunotherapy patients
- Clinical outcome correlation analysis
  
**Data Availability**
- Single-cell data: Available through GEO (accession pending) and dbGaP upon manuscript acceptance
- Bulk RNA-seq data: Publicly available through TCGA and GEO (accession numbers in manuscript)
- Clinical data: De-identified clinical outcomes included with appropriate IRB approvals
  
**Citation**
If you use this code or data, please cite:  
```
@article{mo2025ccl8ccl13,  
  title={CCL8⁺CCL13⁺ tumor-associated macrophages drive early resistance to CAR T cell therapy in large B cell lymphoma},  
  author={Mo, Kelvin C and Kramer, Anne Marijn and others},  
  journal={TBD},  
  year={TBD},  
  publisher={Elsevier}  
}
```
  
**Contact**
- Lead Contact: Dr. Zinaida Good (zinaida@stanford.edu)
- Computational Analysis: Kelvin C. Mo
- Study Oversight: Dr. Zinaida Good, Dr. Crystal L. Mackall, Dr. David B. Miklos
  
**Institutional Affiliations**
- Stanford University School of Medicine
- Center for Cancer Cell Therapy, Stanford Cancer Institute
- Parker Institute for Cancer Immunotherapy
  
**License**
This project is licensed under the MIT License - see the LICENSE file for details
