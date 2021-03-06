# GPCRSampling

## FOLDER DataPreparation
*Contains pymol scripts for data preparation.*
* **FILE alignment_extractor_general_alignment.p**

  * *Pymol script for extracting best rotation points of GPCR helices through alignment of whole GPCR.*
  * Need GPCR's helices borders for executing. Borders could be fond in folder **Data** file **GPCR-TM-table-identity-resis-general.csv**.
* **FILE alignment_extractor_pair_alignment.p**

  * *Pymol script for extracting best rotation points of GPCR helices through pair alignment of GPCR helices. Alinment is processed with exponential smooth with step = 0.5*
  * Need GPCR's helices borders for executing. Borders could be fond in folder **Data** file **GPCR-TM-table-identity-resis-pair.csv**.
* **FILE residues_aligner.p** 

  * *Pymol script which extracts rms_cur data for specific TM and specific GPCR pair.*
  * Need GPCR's helices borders for executing. Borders could be fond in folder **Data** file **GPCR-TM-table-identity-resis-*.csv**.
* **FILE helices-direction-preparation.p**

  * *Pymol script for delices bend direction analysis.*
  * Need GPCR's helices borders for executing. Borders could be fond in folder **Data** file **GPCR-TM-table-directions.csv**.

## FOLDER DataAnalysis
*Contains jupyter notebook and pymol scrips for data analysis.*
* **FILE bend_analysis.ipynb**

  * *Jupyter notebook script for bend points distribution analysis.*
  * Uses **alignment_extractor_general_alignment.p** or **alignment_extractor_pair_alignment.p** from folder **DataPreparation** to obtain data.
  * Uses file **GPCR-TM-table-identity-resis.csv** from folder **Data** as input data.
* **FILE aligned_residues_analysis.ipynb**

  * *Jupyter notebook script for analysis of specific aligneg GPCR pair and specific TM.*
  * Uses **residues_aligner.p** from folder **DataPreparation** to obtain data.
  * Uses file **aligned_residues.csv** from folder **Data** as input data.
* **FILE helices_viewer.p**

  * *Pymol script for specific helix of specific GPCR pair overview.*
  * Need GPCR's helices borders for executing. Borders could be fond in folder **Data** file **GPCR-TM-table-identity-resis-*.csv**.
* **FILE directions_analysis.ipynb**

  * *Jupyter notebook scrip for bend direction analysis*
  * Uses file **helices-direction-preparation.p** from folder **DataPreparation** to obtain data.
  * Uses file **FILE GPCR-TM-table-directions.csv** from folder **Data** as input data
  
## FOLDER Data
*Contains data which is necessary for script execution*
* **FILE GPCR-TM-table-identity-resis-pair.csv**

  * *Input and output of **alignment_extractor_pair_alignment.p**.*
  * *Input of **bend_analysis.ipynb**.*
  * Contains information about best bend parameters obtained through specifit TM alignment.
  * Contain information about GPCR helices borders.
* **FILE GPCR-TM-table-identity-resis-general.csv**

  * *Input and output of **alignment_extractor_general_alignment.p**.*
  * *Input of **bend_analysis.ipynb**.*
  * Contains information about best bend parameters obtained through whole gpcr alignment.
  * Contain information about GPCR helices borders.
* **FILE aligned_residues.csv**

  * *Input for **aligned_residues_analysis.ipynb**.*
  * *Output of **residues_aligner.p**.*
  * Contains information about specific TM alignment.
* **FILE GPCR-TM-table-directions.csv**

  * *Input for **helices-direction-preparation.p**.*
  * *Output of **helices-direction-preparation.p**.*
  * Contains information about helices bend direction.
  
* **GPCR-TM-table-directions-marked-fixed.csv**
  * *Output of **directions_analysis.ipynb**.*
  * Contains information about bend directions with polar coordinates marks.
  * Distances are measured between identical residues.
