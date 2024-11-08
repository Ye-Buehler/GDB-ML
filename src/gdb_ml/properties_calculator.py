from rdkit import Chem
import pandas as pd

class PropertiesCalculator:

    def __init__(self):
        # Initialize any instance variables here, if needed
        pass
    
    # TODO: Compute the complexity score
    def complexity_score(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        divalent_cyclic = 0
        divalent_acyclic = 0
        non_divalent_cyclic = 0
        non_divalent_acyclic = 0

        for atom in mol.GetAtoms():
            degree = atom.GetDegree()
            if degree == 2:
                if atom.IsInRing():
                    divalent_cyclic += 1
                else:
                    divalent_acyclic += 1
            if degree != 2:
                if atom.IsInRing():
                    non_divalent_cyclic += 1
                else:
                    non_divalent_acyclic += 1

        # D1  the fraction divalent cyclic / all cyclic nodes:
        try:
            d1 = round(divalent_cyclic / (divalent_cyclic + non_divalent_cyclic), 2)
        except:
            d1 = 0

        # D2  divalent acyclic /all acyclic nodes:
        try:
            d2 = round(divalent_acyclic / (divalent_acyclic + non_divalent_acyclic), 2)
        except:
            d2 = 0

        # D3  the number of rings:
        number_of_ring = 0
        for s in smiles:
            count_ring_indices = 0
            in_brackets = False
            two_digit_idx = False
            idx = ''
            for char in s:
                if char == '[':
                    in_brackets = True
                elif char == ']':
                    in_brackets = False
                elif char == '%':
                    two_digit_idx = True
                elif not in_brackets and char.isdigit():
                    if two_digit_idx:
                        if len(idx) == 0:
                            idx += char
                        elif len(idx) >= 1:
                            count_ring_indices += 1
                            idx = ''
                            two_digit_idx = False
                    else:
                        count_ring_indices += 1
            ring_part = count_ring_indices / 2
            number_of_ring = number_of_ring + ring_part
            d3 = int(number_of_ring)

        return (d1, d2, d3)


    # TODO: Compute the metric validity and canonicalize the SMILES
    def validity(self, df):

        #to check SMILES validity with handling for NaN and non-string values
        def is_valid_smiles(smiles):
            if isinstance(smiles, str):  # Ensure the entry is a string
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
                    return canonical_smiles
                except: 
                    return False
            return False
        
        df['Is_Valid'] = df['SMILES'].apply(is_valid_smiles)
        false_count = (~df['Is_Valid'].astype(bool)).sum()

        true_percentage = (len(df)-false_count) / len(df) * 100
        df_valid = df[df['Is_Valid'].astype(bool)].reset_index(drop=True)

        return true_percentage, df_valid
    

    # TODO: Compute the metric uniqueness
    def uniqueness(self, df):
        # Get unique rows in the DataFrame
        df_unique = df[['SMILES']].drop_duplicates().reset_index(drop=True)

        uniqueness = df_unique.shape[0]
        return uniqueness, df_unique


    # TODO: Compute the metric novelty
    def novelty(self, df_unique, df_train):
        # Perform an anti-join to get rows in df_unqiue that are not in df_train
        df_novel = pd.merge(df_unique, df_train, how='left', indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)

        novelty = df_novel.shape[0]
        return novelty, df_novel


    # TODO: Compute the metric presence
    def presence(self, df_unique, all_mols_list):
        # Count the occurrences of values in df_unqiue that are in all_mols_list
        presence = df_unique['SMILES'][df_unique['SMILES'].isin(all_mols_list)].reset_index(drop=True)
        presence_count = df_unique['SMILES'].isin(all_mols_list).sum()
        return presence_count, presence
    

    # TODO: To check the ones not in GDB13s but still valid
    def non_presence(self, df_unique, all_mols_list):
        # Find entries in df_unique that are NOT in all_mols_list
        non_presence = df_unique['SMILES'][~df_unique['SMILES'].isin(all_mols_list)].reset_index(drop=True)
        non_presence_count = (~df_unique['SMILES'].isin(all_mols_list)).sum()
        return non_presence_count, non_presence
    

    # TODO: HAC calculation
    def hac(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        size = mol.GetNumHeavyAtoms()

        return size