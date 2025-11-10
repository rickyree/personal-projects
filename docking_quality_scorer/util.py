import pandas as pd
import numpy as np
from math import *
import os 
from Bio.PDB import PDBParser 
import anarci 
from difflib import SequenceMatcher 
import matplotlib.pyplot as plt 
from Bio.PDB import NeighborSearch 
from collections import Counter 



aa_keys = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}

#create a dataframe for the csv file (raw) containing the heatmap 
def create_heatmap(raw, by_CDR = False): 
    df = pd.read_csv(raw) 
    if by_CDR: 
        region = [x for x in df['region'].unique()]
    residues = [x for x in df['abres'].unique()]
    df.dropna() 
    #dict containing all heatmaps for each CDR/FR region: 
    CDR_heatmaps = {} 

    #if the heatmap is not assorted by CDRS: 
    if not by_CDR: 
        heatmap = pd.DataFrame(columns = df['abres'].unique())
        for x,y in enumerate(residues):
            heatmap.loc[x] = [0.0 for x in residues]
        heatmap.index = residues

        for index, row in df.iterrows():
            heatmap.loc[row['abres'], row['agres']] = heatmap.loc[row['abres'], row['agres']] + row['score']

        heatmap = round(heatmap, 1).astype(float) 
        return heatmap 

    #if the heatmap is assorted by CDRs: 
    elif by_CDR: 
        #function for creating a heatmap per CDR: 
        def create_CDR_heatmap(CDR_label): 
            CDR_heatmap = pd.DataFrame(columns = df['abres'].unique())
            for x,y in enumerate(residues):
                CDR_heatmap.loc[x] = [0.0 for x in residues]
            CDR_heatmap.index = residues 

            for index, row in df[df['region'] == CDR_label].iterrows():
                CDR_heatmap.loc[row['abres'], row['agres']] = CDR_heatmap.loc[row['abres'], row['agres']] + row['score'] 
            
            CDR_heatmap = CDR_heatmap.fillna(0)
            CDR_heatmap = round(CDR_heatmap, 1).astype(float) 
            return CDR_heatmap 
        
        #create a heatmap for each CDR and add to dict 
        CDR_heatmaps['CDRL1'] = create_CDR_heatmap('CDR-L1')
        CDR_heatmaps['CDRL2'] = create_CDR_heatmap('CDR-L2')
        CDR_heatmaps['CDRL3'] = create_CDR_heatmap('CDR-L3')
        CDR_heatmaps['CDRH1'] = create_CDR_heatmap('CDR-H1')
        CDR_heatmaps['CDRH2'] = create_CDR_heatmap('CDR-H2') 
        CDR_heatmaps['CDRH3'] = create_CDR_heatmap('CDR-H3') 
        CDR_heatmaps['None'] = create_CDR_heatmap('') 

        return CDR_heatmaps 


#search path for output file (should be formatted as "output.csv") and assigning docking scores to each pdb: 
def searchdir_output(mypath, docking_scores): 
    for (root, _, files) in os.walk(mypath): 
        for file in files: 
            if file == 'output.csv': 
                for index, row in pd.read_csv(os.path.join(root, file)).iterrows(): 
                    if row['model'][-4:] == '.pdb': 
                    #and row['model'] not in ['processed_NTF3_NTF3-FAB_AZ13490692_5fear_1.pdb', 'processed_FAb-HEL_Lyso0025-Lyso0083_None_0beth_1.pdb']: 
                    #if len(row['pdb']) == 4: 
                        docking_scores[row['model']] = row['DockQ'] 
                        #docking_scores[row['pdb'] + '.pdb'] = row['-log10(kd)'] 
    
    return docking_scores


#search path for pdb files and make an array of the files for analysis: 
def searchdir_files(mypath): 
    filelist = [] 
    for (root, _, files) in os.walk(mypath): 
        for file in files: 
            if (len(file) >= 4) and (file[-3:] == 'pdb'): 
                    filelist.append(os.path.join(root, file)) 
    
    return filelist 



#classifying PDBs into heavy, light and antigen chains and identifying CDR regions: 
def pdb_label(file): 
    #parse pdb file 
    pdbparser = PDBParser(QUIET = True) 
    structure = pdbparser.get_structure('structure1', file) 

    '''if len([chain for chain in structure.get_chains()]) < 3: 
        print("Minimum of 3 chains required for heavy, light and antigen chain. ") 
        return '', '', '', '', '', True  '''

    #concatenated sequence of all residues in the structure: 
    seq = '' 

    #analyzing sequence with ANARCI:    
    for residue in structure.get_residues(): 
        if residue.get_resname() in list(aa_keys.keys()): 
                seq = ''.join([seq, aa_keys[residue.get_resname()]])
    sequences = [('asdf', seq)]

    numbering, alignment_details, _ = anarci.anarci(sequences, scheme = 'IMGT', output = False) 
    return structure, numbering, alignment_details 


def chain_and_CDR_allocation(structure, numbering, alignment_details): 
    #df containing information of all the residues, their chains (H/L), position, IMGT number and CDR region: 
    residues_sorted = pd.DataFrame(columns = ['residue', 'chain', 'pos', 'IMGT', 'CDR']) 

    #allocate residue and its position 
    for residue in structure.get_residues(): 
        if residue.get_resname() in list(aa_keys.keys()): 
            residues_sorted.loc[len(residues_sorted), 'pos'] = residue.get_id()[1] 
            residues_sorted.loc[len(residues_sorted)-1, 'residue'] = aa_keys[residue.get_resname()] 

    residue = [n for n in residues_sorted['residue']] 

    #to be concatenated with H or L sequence: 
    H = '' 
    L = ''

    #for each chain (H and L): 
    for x in range(len(numbering[0])): 
        #store IMGT numbering: 
        IMGTnumbers = [n[0][0] for n in numbering[0][x][0]] 
        #store residue if applicable to IMGT number: 
        IMGTresidues = [n[1] for n in numbering[0][x][0]] 
        
        if H == '': 
            hflag = True 
        else: 
            hflag = False 
        
        if L == '': 
            lflag = True 
        else: 
            lflag = False 
        
        #set r as starting residue number of H or L chain: 
        #set a and n as 0 to work its way through the IMGT labelings: 
        r = alignment_details[0][x]['query_start'] 
        n = 0
        a = 0 

        #until r reaches end of chain: 
        while r < alignment_details[0][x]['query_end']: 
            if IMGTresidues[a] == residue[r]: 

                #add IMGT and chain label to residues_sorted table: 
                residues_sorted.loc[r, 'IMGT'] = IMGTnumbers[n] 
                residues_sorted.loc[r, 'chain'] = alignment_details[0][x]['chain_type'] 

                #if K change to L 
                if residues_sorted.loc[r, 'chain'] == 'K': 
                    residues_sorted.loc[r, 'chain'] = 'L'

                #concatenate the residue with the H/L sequence: 
                if residues_sorted.loc[r, 'chain'] == 'H' and hflag == True: 
                    H = ''.join([H, residues_sorted.loc[r, 'residue']])
                if residues_sorted.loc[r, 'chain'] == 'L' and lflag == True: 
                    L = ''.join([L, residues_sorted.loc[r, 'residue']]) 
                                
                #according to IMGT numbering: 
                if 1 <= IMGTnumbers[n] <= 26: 
                    residues_sorted.loc[r, 'CDR'] = 'FR1' 
                elif 27 <= IMGTnumbers[n] <= 38: 
                    residues_sorted.loc[r, 'CDR'] = 'CDR1' 
                elif 39 <= IMGTnumbers[n] <= 55: 
                    residues_sorted.loc[r, 'CDR'] = 'FR2' 
                elif 56 <= IMGTnumbers[n] <= 65: 
                    residues_sorted.loc[r, 'CDR'] = 'CDR2' 
                elif 66 <= IMGTnumbers[n] <= 104: 
                    residues_sorted.loc[r, 'CDR'] = 'FR3' 
                elif 105 <= IMGTnumbers[n] <= 117: 
                    residues_sorted.loc[r, 'CDR'] = 'CDR3' 
                elif 118 <= IMGTnumbers[n] <= 128: 
                    residues_sorted.loc[r, 'CDR'] = 'FR4' 


                #move on to next residue: 
                r += 1 
                a += 1 
                n += 1 
        
            else: 
                a += 1 
                n += 1 
    return H, L, residues_sorted


def label_chains_and_CDR(structure, H, L, residues_sorted): 
    #labelling each chain in PDB file as antibody or antigen: 

    #dict containing which chain corresponds to H or L. 
    #in the format {H: A, L: B} if the chains in PDB are labelled as A and B
    antibody_chains = {'H' : [], 'L' : []} 
    #array containing all the chains that belong to the antigen 
    antigen_chains = [] 
    #dict containing identified CDR residues 
    CDR_residues = {}
    #warning if there is no H or L chain: 
    warning = False 

    #store residue positions of each CDR: 
    def label_CDR_residues(chain, CDR_residues): 
        for n in range(1, 4): 
            CDR_residues[f'CDR{chain}{n}'] = [x for x in residues_sorted['pos'][residues_sorted.index[(residues_sorted['chain'] == chain) & (residues_sorted['CDR'] == f'CDR{n}')].tolist()]] 
    
    def label_FR_residues(chain, CDR_residues): 
        for n in range(1, 4): 
            CDR_residues[f'FR{chain}{n}'] = [x for x in residues_sorted['pos'][residues_sorted.index[(residues_sorted['chain'] == chain) & (residues_sorted['CDR'] == f'FR{n}')].tolist()]] 

    #for each chain, if the chain matches identified H/L chain, label it as such and add to antibody_chains, otherwise add to antigen_chains: 
    #print('H: ', H) 
    #print('L: ', L)
    for chain in structure.get_chains(): 
        seq = '' 
        for residue in chain.get_residues(): 
            if residue.get_resname() in list(aa_keys.keys()): 
                seq = ''.join([seq, aa_keys[residue.get_resname()]]) 

        if len(seq) > 0: 
            if (SequenceMatcher(None, seq, L[:len(seq)]).ratio() > 0.95) or (L[3:-3] in seq): 
                antibody_chains['L'].append(chain.get_id()) 
                label_CDR_residues('L', CDR_residues) 
                label_FR_residues('L', CDR_residues) 
                #print('added into L: ', chain.id, seq)

            elif (SequenceMatcher(None, seq, H[:len(seq)]).ratio() > 0.95) or (H[3:-3] in seq): 
                antibody_chains['H'].append(chain.get_id()) 
                label_CDR_residues('H', CDR_residues) 
                label_FR_residues('H', CDR_residues) 
                #print('added into H:', chain.id, seq)

            else: 
                antigen_chains.append(chain.get_id()) 
                #print('added into antigen: ', chain.id, seq)

    #return warning if there are no H or L chain and move on to next PDB file: 
    if ([x for x in antibody_chains.keys()].count('H') == 0) or ([x for x in antibody_chains.keys()].count('L') == 0): 
        print('Heavy/light chains not correctly assigned. ') 
        warning = True 
    return structure, antibody_chains, antigen_chains, CDR_residues, residues_sorted, warning


#obtaining contacting residue pairs: 
def get_contact_pairs(CDR_residues, antibody_chains, antigen_chains, structure, distance, only_H3): 
    #dataframe summarizing all the residue pairs, their chain, CDR region, complementarity scores, contacting atoms, coordinates and distance: 
    contact_pairs = [] 
    hindrance_counts = 0 

    atoms_to_exclude = ['H', 'OH', 'NH1', 'NH2', 'CH']

    #identify contact pairs using NeighborSearch per each antibody and antigen chain combination: 
    def get_contact_pairs_by_chain(heavy_or_light, antibody_chain, antigen_chain, n, contact_pairs): 
        
        #residue, atom 1: antibody chain 
        #residue, atom 2: antigen chain 

        #store all non-H atoms in the antigen chain under atom2s array: 
        atom2s = [] 
        for residue2 in structure[0][antigen_chain].get_residues(): 
            if residue2.get_resname() in list(aa_keys.keys()): 
                for atom2 in residue2.get_atoms(): 
                    #exclude hydrogens 
                    if not atom2.get_id().startswith('H') and not (atom2.get_id()[0].isdigit() and atom2.get_id()[1] == 'H'): 
                        atom2s.append(atom2) 
        
        #create NeighborSearch class for the antigen chain: 
        search_pairs = NeighborSearch(atom2s) 

        #store all residues in the antibody chain under residue1s array as tuple (residue1, CDR region): 
        
        residue1s = [] 
        for residue in structure[0][antibody_chain].get_residues(): 
            if only_H3: 
                if residue.get_id()[1] in CDR_residues['CDRH3']: 
                    residue1s.append((residue, 'CDRH3')) 
            elif not only_H3: 
                if (residue.get_id()[1] in CDR_residues[f'CDR{heavy_or_light}{n}']): 
                    residue1s.append((residue, f'CDR{heavy_or_light}{n}'))
        
        #iterate through each atom1 in residue1 to obtain contacting atom2s: 
        for residue1, CDR in residue1s: 
            residues_found = [] 
            for atom1 in residue1.get_atoms(): 
                #exclude hydrogens 
                if not atom1.get_id().startswith('H') and not (atom1.get_id()[0].isdigit() and atom1.get_id()[1] == 'H'): 
                    temp = search_pairs.search(center = atom1.get_coord(), radius = distance, level = 'R') 
                    if temp != []: 
                        residues_found.extend(temp) 
            residues_found = set(residues_found) 

            #add each residue1 and residue2 to contact_pairs: 
            for residue_found in residues_found: 
                contact_pairs.append((''.join((residue1.get_resname()[:3], residue_found.get_resname()[:3])), CDR))    
        
        return contact_pairs 
    
    #identify pairs and add to contact_pairs for each antibody chain and CDR 1 to 3
    for antigen_chain in antigen_chains: 
        
        for antibody_chain in antibody_chains['H']: 
            for n in range(1,4): 
                contact_pairs = get_contact_pairs_by_chain('H', antibody_chain, antigen_chain, n, contact_pairs)
        for antibody_chain in antibody_chains['L']: 
            for n in range(1,4): 
                contact_pairs = get_contact_pairs_by_chain('L', antibody_chain, antigen_chain, n, contact_pairs) 
    return contact_pairs, hindrance_counts



#scoring residue pairs: 
def residue_score(CDR_residues, contact_pairs, heatmap, by_CDR, only_H3): 
    x = 0
    score = 0 
    label = 'incorrect'

    #score_by_CDR stores information for each CDR region as tuple (total score, frequency of residues) 
    score_by_CDR = {} 

    for cdr in ['CDRH1', 'CDRH2', 'CDRH3', 'CDRL1', 'CDRL2', 'CDRL3']: 
        score_by_CDR[cdr] = [0, 0]

    #for each pair in contact pairs: 
    while x < len(contact_pairs): 
        #fetch name of each residue of the xth pair: 
        residue1 = contact_pairs[x][0][:3]
        residue2 = contact_pairs[x][0][3:]

        #identify scores for the residues based on the heatmap: 
        if not by_CDR: 
            n = heatmap.loc[residue1,residue2] 
        elif by_CDR: 
            #identify residues in heatmap: 
            n = heatmap[contact_pairs[x][1]].loc[residue1,residue2] 
            if not only_H3: 
                score_by_CDR[contact_pairs[x][1]][0] += n
                score_by_CDR[contact_pairs[x][1]][1] += 1
        
        if not np.isnan(n): 
            #also add to the "score" column in contact_pairs: 
            #contact_pairs.loc[x, 'score'] = n 
            score += round(n, 1) 
        
        #move on to the next contact pair: 
        x+=1 
    
    if score < -4: 
        label = 'incorrect' 
    elif score < 10: 
        label = 'carry forward' 
    else: 
        label = 'carry forward with priority'
    
    return score, contact_pairs, score_by_CDR, label


def export_files(obtained_scores): 
    results = pd.DataFrame(obtained_scores, columns = ['file', 'c-score', 'label']) 
    os.makedirs('results', exist_ok = True)
    results.to_csv('results/results.csv') 


def labels_pichart(obtained_scores): 
    label_counts = Counter([x[2] for x in obtained_scores])

    # Plot pie chart
    plt.figure(figsize=(6, 6))
    plt.pie(
        label_counts.values(),
        labels=label_counts.keys(),
        autopct='%1.1f%%',
        startangle=140,
        colors=['#d73027', '#1a9850', '#91bfdb']  # Optional: assign custom colors
    )
    plt.title('Label Distribution')
    plt.tight_layout()
    plt.savefig('results/pie.png') 
