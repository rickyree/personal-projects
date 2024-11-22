## Welcome! 
This is a scorer for antibody-antigen complexes (.pdb) by assessing each contacting residue pairs with a residue by residue heatmap taken from a publication. To get started, run the file affinity_scorer.py, with directories and variables under the 'variables to modify' section adjusted accordingly. 

## Files  

### affinity_scorer.py 
Takes a folder containing the PDB files and output file (should be named as 'output.csv'). Can be taken from SabDab. Calculates complementarity scores from contacting residue pairs and returns a graph assessing the scores with kd values in the results folder. 

### results 
Folder containing graphs correlating obtained complementarity scores with kd/docking scores. Arranged by publications that published the residue heatmap. Results of running affinity_scorer.py will be stored here. 

### example_dataset 
Sample database containing pdb files and their corresponding kd values. Set this as rootdir in affinity_scorer.py. 

### residuemaps 
Folder containing residue heatmaps for each publication. Set one of the csv files inside as heatmap_csv in affinity_scorer.py. 
* residuemap_akbar.csv: arranged by CDR 
* residuemap_li.csv: not arranged by CDR 
* residuemap_ramaraj: not arranged by CDR 

### miscellaneous 
A bunch of other files for practice. 

## Citations for heatmaps 
* Akbar, R. et al. (2021) ‘A compact vocabulary of paratope-epitope interactions enables predictability of antibody-antigen binding’, Cell Reports, 34(11), p. 108856. doi:10.1016/j.celrep.2021.108856. 
* Ramaraj, T. et al. (2012) ‘Antigen–antibody interface properties: Composition, Residue Interactions, and features of 53 non-redundant structures’, Biochimica et Biophysica Acta (BBA) - Proteins and Proteomics, 1824(3), pp. 520–532. doi:10.1016/j.bbapap.2011.12.007. 
* Li, J. et al. (2024) ‘Development and experimental validation of computational methods for human antibody affinity enhancement’, Briefings in Bioinformatics, 25(6). doi:10.1093/bib/bbae488. 




