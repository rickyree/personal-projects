import pandas as pd
from math import *
import os 
import argparse 
from util import *




"""
---Variables to modify---

rootdir: directory containing the PDB files and output file (should be named as "output.csv"). Can be taken from SabDab. 

heatmap_csv: csv file containing residue by residue heatmap. 
- columns should be labeled as "abres" for antibody residues, "agres" for antigen residues, "region" for CDR region, "score" for affinity score 
- residues should be in 3 letter format e.g. ASP 
- CDR regions should be in the format e.g. CDR-H1 or HFR1 

result_name: name to save the results file in. Saved under the "results" folder. 

only_H3: set to True if the analysis is only to be conducted on H3 regions. 
"""


def main(): 
    print('Welcome! \nLoading...') 

    aa_keys = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}

    #parser for command line interface 
    parser = argparse.ArgumentParser(description = 'Evaluate docking complexes by c-scores. ') 

    required_args = parser.add_argument_group('required') 

    required_args.add_argument('-b', '--base', type = str, required = True, help = 'Path where PDB files are stored.\n') 

    parser.add_argument('-hm', '--heatmap_csv', type = str, default = 'residuemap/residuemap.csv', help = 'Path where residue heatmap is stored. ') 

    parser.add_argument('-d', '--distance', type = int, default = 5, help = 'Distance cutoff of contact pairs in Angstroms (default = 5A). ') 

    parser.add_argument('-o', '--only_H3', action = 'store_true', help = 'Set to true if the analysis is only to be conducted on H3 regions (default = false). ') 

    #parser.add_argument('--results_name', type = str, help = 'What to name the results file. End with file format e.g. .png') 


    args = parser.parse_args() 

    base = args.base 
    heatmap_csv = args.heatmap_csv 
    distance = args.distance 
    only_H3 = args.only_H3 
    results_name = 'results.png'


    by_CDR = True
            
    #heatmap: df storing all the residue-residue complementarity scores 
    #obtained_scores: array to store the complementarity and DockQ scores 
    #root: base directory of pdb files 
    #scores_all: df to store scores by CDR 
    #filelist: list of PDB files to analyse
    heatmap = create_heatmap(heatmap_csv, by_CDR)
    obtained_scores = []
    blanks = []
    filelist = [] 
    scores_all = pd.DataFrame(columns = ['CDRH1', 'CDRH2', 'CDRH3', 'CDRL1', 'CDRL2', 'CDRL3', 'DockQ']) 
        

     
    filelist = searchdir_files(base) 

    for file in filelist: 
        print(os.path.basename(file))
        

        #classifying PDBs into heavy, light and antigen chains and identifying CDR regions: 
        structure, numbering, alignment_details = pdb_label(file) 
        
        #classifying residues into CDRs and IMGT numbering: 
        #H, L: sequence of H and L chains 
        #residues_sorted: df containing information of all the residues, their chains (H/L), position, IMGT number and CDR region: 
        H, L, residues_sorted = chain_and_CDR_allocation(structure, numbering, alignment_details) 

        #labeling chain and CDR regions and skip to next file if no H or L chain: 
        structure, antibody_chains, antigen_chains, CDR_residues, residues_sorted, warning = label_chains_and_CDR(structure, H, L, residues_sorted) 
        if warning: 
            print(antibody_chains, CDR_residues) 
            continue 
        
        #obtaining contacting residue pairs: 
        contact_pairs, hindrance_counts = get_contact_pairs(CDR_residues, antibody_chains, antigen_chains, structure, distance, only_H3) 
        if len(contact_pairs) < 1: 
            print('no contact pairs found! ') 
            continue 
        
        #scoring residue pairs: 
        score, contact_pairs, score_by_CDR, label = residue_score(CDR_residues, contact_pairs, heatmap, by_CDR, only_H3) 
        
        
        if not only_H3 and by_CDR: 
            count = 0
            for x in score_by_CDR.keys(): 
                if score_by_CDR[x][1] >= 1: 
                    count += 1
            if count >= 1 and -100 < score < 100: 
                #saves scores in a tuple as such: (complementarity score, dockQ score, name of protein structure) 
                score = round(float(score), 2)
                obtained_scores.append((os.path.basename(file), score, label)) 
                print('c-score: ', score, '\n') 

        elif only_H3 or not by_CDR: 
            #saves scores in a tuple as such: (complementarity score, dockQ score, name of protein structure) 
            score = round(float(score), 2)
            obtained_scores.append((os.path.basename(file), score, label)) 
            print('c-score: ', score, '\n') 
    
    #export scores and labels as csv: 
    export_files(obtained_scores)

    #export pi chart of labels: 
    labels_pichart(obtained_scores) 

    


if __name__ == '__main__': 
    main()  




