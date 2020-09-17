# Skyline-Protein-Turnover
Protein Turnover tool in Skyline

The Skyline Protein Turnover tool uses multiple R scripts to calculate the rates of protein turnover from mass spectrometry proteomics experiments.

Step0
Input File: Skyline report csv file (e.g. 2020_0723_OCR_turnover_v1.csv) 
Input Parameters: diet.enrichment, min.abundance, resolution, p.tolerance, label (in progress- currently the only label allowed is leucine)
Output Files: 
1) multiple leucine peptides data file (e.g. Step0_Data_Output_Skyline_multileucine_peps_test.csv)
2) single leucine peptide data file (e.g.Step0_Data_Output_Skyline_singleleucine_peps_test.csv)
This script is the main workhorse as it corrects for the naturally occurring heavy isotopes of Hydrogen, Carbon, Oxygen, Nitrogen and Sulfur. This is the most compute-intensive step in this tool.


Step1
Step2

Step3
Input Files:
1) multiple leucine peptides data file (e.g. Step0_Data_Output_Skyline_multileucine_peps_test.csv)
2) single leucine peptide data file (e.g.Step0_Data_Output_Skyline_singleleucine_peps_test.csv)

Step4
Step5
