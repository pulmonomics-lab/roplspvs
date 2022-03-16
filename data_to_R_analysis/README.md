# Description of datafiles

Data was downloaded from Metabolights https://www.ebi.ac.uk/metabolights/MTBLS136/files. 
"m_metabolite_profiling_mass_spectrometry_v2_maf.tsv" containing datamatrix was renamed into "s_MTBLS136_datamatrix.txt" and edited as follows:
1. Metabolite descriptions were removed.
2. Matrix was transposed to have metabolite names as column names in first row.
3. Samples that did not have a group identification was removed.
4. Special characters in variable names were replaced by underscore or text.

"s_MTBLS136.txt" containing metadata was renamed into "sampleID.txt" and edited as follows:
1. special characters were replaced by underscore
2. "Sample Name" column was moved to column 1
3. Samples were ordered to correspond to datamatrix
4. "<=" was replaced with "less or equal to" in secondary ID description
5. Samples that did not have a group identification was removed
