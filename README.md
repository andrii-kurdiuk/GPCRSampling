# GPCRSampling
## FOLDER DataPreparation
*Contains pymol scripts for data preparation*
* **FILE alignment_extractor.p**
...Pymol script for extracting best rotation points of GPCR helices
** need GPCR's helices borders for executing
## FILE GPCR-TM-table-identity-resis.csv
*input and output of **alignment_extractor.p** - file with GPCR structural information*
* contain information about best rotation parameters
* contain information about GPCR helices borders
## FILE residues_aligner.p
*extract rms_cur data for specific TM and specific GPCR pair*
* need GPCR's helices borders for executing
## FILE bend_analysis.ipynb
*analysis of bend point distribution*
* use **alignment_extractor.p** to obtain data
## FILE aligned_residues_analysis.ipynb
*File for analysis of specific aligneg GPCR pair and specific TM*
* use **residues_aligner.p** to obtain data


