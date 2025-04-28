## Welcome! 
This is a scoring system for assessing antibody-antigen docking quality by assessing the contact preferences (c-score) of each contacting residue pairs. To get started, run the file cscore_generator.py, specifying the `-b` option with the path to the PDB files. 

## Files  

### cscore_generator.py 
Takes a folder containing the PDB files of docking complexes to generate c-scores for. Calculates c-scores from contacting residue pairs and returns labels for each complex to discard or carry forward to further analysis in the results folder. 

### util.py 
File containing all internal functions 

### results 
Results of running cscore_generator.py will be stored here. After running, 2 files are generated: 
* generated_scores.csv: csv file of all c-scores and labels for each complex. 
* pie.png: pie chart showing the proportions of complex labels.  

### example_files 
Sample database containing PDB files to generate c-scores for. Set this with the `-b` option in cscore_generator.py. Note these are SAbDab structures and not docking models. 

### residuemap.csv 
CSV file containing contact preference residue heatmaps. 

### miscellaneous 
A bunch of other files for practice. 

## Citation for contact preference heatmap 
* Akbar, R. et al. (2021) ‘A compact vocabulary of paratope-epitope interactions enables predictability of antibody-antigen binding’, Cell Reports, 34(11), p. 108856. doi:10.1016/j.celrep.2021.108856. 



