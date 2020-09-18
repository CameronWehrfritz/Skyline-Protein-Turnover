# Skyline-Protein-Turnover
Protein Turnover tool in Skyline

The Skyline Protein Turnover tool uses multiple R scripts to calculate the rates of protein turnover from mass spectrometry proteomics experiments.

Step1
Description: This script is the main workhorse as it corrects for the naturally occurring heavy isotopes of Hydrogen, Carbon, Oxygen, Nitrogen and Sulfur. This is the most compute-intensive step in this tool.
Input File: Skyline report csv file (e.g. 2020_0723_OCR_turnover_v1.csv) 
Input Parameters: diet.enrichment, min.abundance, resolution, p.tolerance, label (in progress- currently the only label allowed is leucine), isotope dot product filter, average turnover score filter
Output Files: 
1) multiple leucine peptides data file (e.g. Step0_Data_Output_Skyline_multileucine_peps_test.csv)
2) single leucine peptide data file (e.g.Step0_Data_Output_Skyline_singleleucine_peps_test.csv)

Step2
Description: Optional step to calculate the x-intercepts of the turnover regressions, using the percentage newly synthesized values of each protein, calculated step 1.
Input File:
Input Parameters:
Output Files: 

Step3
Description: Performs the regressions of percentage newly synthesized proteins (calculated in step 1). The x-intercepts of the regressions are fixed to values calculated in Step 2. If step 2 was note used, then regressions go through the origin. Generates a report of the slope values and statistical summaries of the regressions.
Input File: 
Input Parameters: x-intercept values (read from step 2, or defaults to 0)
Output Files: 1) CSV report containing slope values, half-lives, x-intercepts, errors, p-values associated with each regression, 2) PDF file containing regressions for every protein and condition

Step4
Description: Performs pairwise statistical comparisons of all the regressions between treatment groups and generates statistical reports.
Input Files: CSV file outputted from step 3.
1) multiple leucine peptides data file (e.g. Step0_Data_Output_Skyline_multileucine_peps_test.csv)
2) single leucine peptide data file (e.g.Step0_Data_Output_Skyline_singleleucine_peps_test.csv)
Input Parameters: confidence interval filters (not implemented yet)
Output Files: 1) CSV report containing the half-lives, log fold changes between treatment groups, p-values, q-values and other statistical summary information


