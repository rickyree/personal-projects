## Welcome! 
This is a scoring system for assessing antibody-antigen docking quality by assessing the contact preferences (c-score) of each contacting residue pairs. Contact preferences are based on interaction residue frequencies across Ab-Ag complexes in the AbDb database (http://www.abybank.org/abdb/). 

To get started, create a virtual environment with the required packages in requirements.txt installed. Note: this doesn't include anarci, which needs to be separately installed as per instructions in the repository (https://github.com/oxpig/ANARCI). **You can run the command `source create_environment.sh` which does all of this for you and creates an environment with anarci.** 

Then, run the file cscore_generator.py, specifying the `-b` option with the path to the PDB files (e.g. `python cscore_generator.py -b path_with_example_files`). `-h` for details. 

The c-score for each Ab-Ag complex and their accuracy label (discard, carry forward, or carry forward with priority) will be available as a csv file in the results folder. 


## Files  

### create_environment.sh
Creates and enters into a virtual environment with all dependent packages installed. 

### cscore_generator.py 
Takes a folder containing the PDB files of docking complexes to generate c-scores for. Calculates c-scores from contacting residue pairs and returns labels for each complex to discard or carry forward to further analysis in the results folder. 

### miscellaneous 
A bunch of other files for practice. 

### results 
Results of running cscore_generator.py will be stored here. After running, 2 files are generated: 
* generated_scores.csv: csv file of all c-scores and labels for each complex. 
* pie.png: pie chart showing the proportions of complex labels.  

### path_with_example_files.zip 
Sample database containing PDB files to generate c-scores for. Set this with the `-b` option in cscore_generator.py. Note these are SAbDab structures and not docking models. 

### residuemap
folder with CSV file containing contact preference residue heatmaps. 

### util.py 
File containing all internal functions. 





