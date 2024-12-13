
import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import RDConfig
from rdkit.Chem import QED
from rdkit.Chem.Descriptors import qed
import openbabel as ob
from openbabel import pybel
import sascorer
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))


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

        df_valid = pd.DataFrame({'SMILES': df_valid['Is_Valid'], 'Log Probs': df_valid['Log Probs']})
        
        df_valid = df_valid.dropna(subset=['SMILES']).reset_index(drop=True)

        return true_percentage, df_valid
    

    # TODO: Compute the metric uniqueness
    def uniqueness(self, df):
        # Get unique rows in the DataFrame
        df_unique = df.drop_duplicates(subset=['SMILES']).reset_index(drop=True)

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
        presence = df_unique[df_unique['SMILES'].isin(all_mols_list)].reset_index(drop=True)
        presence_count = df_unique['SMILES'].isin(all_mols_list).sum()
        return presence_count, presence
    

    # TODO: To check the ones not in GDB13s but still valid
    def non_presence(self, df_unique, all_mols_list):
        # Find entries in df_unique that are NOT in all_mols_list
        non_presence = df_unique[~df_unique['SMILES'].isin(all_mols_list)].reset_index(drop=True)
        non_presence_count = (~df_unique['SMILES'].isin(all_mols_list)).sum()
        return non_presence_count, non_presence
    

    # TODO: HAC calculation
    def hac(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        size = mol.GetNumHeavyAtoms()

        return size
    

    # TODO: AlogP calculation
    def alogp(self, smiles):
        molecule = Chem.MolFromSmiles(smiles)
        alogp = Descriptors.MolLogP(molecule)

        return alogp
    

    # TODO: Calculate the FDV
    def divalent_nodes_fraction(self, smiles):
        mol=Chem.MolFromSmiles(smiles)
        atom_number = 0
        divalent_node = 0
        
        for atom in mol.GetAtoms():
            atom_number += 1
            degree = atom.GetDegree()

            if degree == 2:
                divalent_node += 1 
            else:
                continue

        divalent_ratio = round(divalent_node/atom_number,2)
        return divalent_ratio
    

    # TODO: Calculate the SAS
    def sascore(self, smiles):
        mol=Chem.MolFromSmiles(smiles)
        sas_score = sascorer.calculateScore(mol)
        return sas_score


    # TODO: filtration - max 3 rings
    def threeringcheck(self, smiles):
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
            
        if number_of_ring <= 3:
            return True
        else:
            return False


    # TODO: filtration - no atom in 3 rings
    def has_atom_in_three_rings(self, smiles):
        mol = pybel.readstring("smi", smiles)

        # Perform ring perception
        mol.OBMol.AddHydrogens()
        mol.OBMol.PerceiveBondOrders()

        # Get the number of rings and their sizes
        ring_list = []
        for ring in mol.OBMol.GetSSSR():
            #print(ring)
            ring_list.append(ring)

        # Print the ring sizes
        #print("Ring list:", ring_list)

        atom_ring_count = {}
        # Iterate over the rings
        for ring in ring_list:
            # Iterate over the atoms in the ring
            for atom_idx in ring._path:
                if atom_idx in atom_ring_count:
                    atom_ring_count[atom_idx] += 1
                else:
                    atom_ring_count[atom_idx] = 1
        atom_ring_count

        for atom, ring_count in atom_ring_count.items():
            if ring_count >= 3:
                return True

        return False


    # TODO: filtration - has 3 or 4 membered rings
    def has_small_rings(self, smiles):
        mol = pybel.readstring("smi", smiles)

        # Perform ring perception
        mol.OBMol.AddHydrogens()
        mol.OBMol.PerceiveBondOrders()

        # Get the number of rings and their sizes
        small_ring_list = []
        for ring in mol.OBMol.GetSSSR():
            # Check if the ring has not more than 4 atoms (small ring)
            if len(ring._path) <= 4:
                small_ring_list.append(ring)
        
        # Check if the list is empty
        if not small_ring_list:
            #print("The small_ring_list is empty.")
            return False

        else:
            #print("The small_ring_list  is not empty.")
            return True
    
    
    # TODO: filtration - has NO or NN in non-aromatic rings or acyclic bonds
    def if_NO_NN_in_non_aromatic_ring_or_acyclic(self, smiles):

        def isRingAromatic(m, bondRing):
            for id in bondRing:
                if not m.GetBondWithIdx(id).GetIsAromatic():
                    return False
            return True
        
        mol = Chem.MolFromSmiles(smiles)

        for bond in mol.GetBonds():

            atom1_idx = bond.GetBeginAtomIdx()
            atom2_idx = bond.GetEndAtomIdx()

            atom1 = mol.GetAtomWithIdx(atom1_idx)
            atom2 = mol.GetAtomWithIdx(atom2_idx)
            symbol1 = atom1.GetSymbol()
            symbol2 = atom2.GetSymbol()
        
            if bond.IsInRing() == True:
                
                if atom1.GetIsAromatic() == False or atom2.GetIsAromatic() == False:
                    #print(symbol1, symbol2)
                    
                    if symbol1 == 'N' and symbol2 == 'O':
                        #print ('find a NO', symbol1, symbol2)
                        return True

                    elif symbol1 == 'O' and symbol2 == 'N':
                        #print ('find a ON', symbol1, symbol2)
                        return True

                    elif symbol1 == 'N' and symbol2 == 'N':
                        #print ('find a NN', symbol1, symbol2)
                        return True
                    
            elif bond.IsInRing() == False:

                if symbol1 == 'N' and symbol2 == 'O':
                    #print ('find a NO', symbol1, symbol2)
                    return True

                elif symbol1 == 'O' and symbol2 == 'N':
                    #print ('find a ON', symbol1, symbol2)
                    return True

                elif symbol1 == 'N' and symbol2 == 'N':
                    #print ('find a NN', symbol1, symbol2)
                    return True
        return False
    

    # TODO: filtration - has N in three rings
    def if_contain_N3ring(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        ringinfo = mol.GetRingInfo()
        for ring in ringinfo.AtomRings():
            if len(ring) == 3:
                #print(ring)
                idx1 = ring[0]
                idx2 = ring[1]
                idx3 = ring[2]
                #print(idx1,idx2,idx3)

                atom1=mol.GetAtomWithIdx(idx1).GetSymbol()
                atom2=mol.GetAtomWithIdx(idx2).GetSymbol()
                atom3=mol.GetAtomWithIdx(idx3).GetSymbol()
                three_member_ring = atom1, atom2, atom3
                #print(three_member_ring)
                if 'N' in three_member_ring:
                    return True
        return False
    

    # TODO: filtration - has N in three rings
    def if_contain_OCO(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        patt = Chem.MolFromSmarts('OCO')
        try:
            hit_ats = list(mol.GetSubstructMatch(patt))
            hit_bonds = []
            for bond in patt.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())
            return True
        except:
            return False
        

    # TODO: filtration - non-aromatic double bond filter
    def non_aromatic_double_bond_filter(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        for bond in mol.GetBonds():
            bondtype = bond.GetBondType()
            if bondtype == 2: # aromatic bond cannot be 2
                #print("It is a double bond")
                mol_begin=bond.GetBeginAtomIdx()
                mol_end=bond.GetEndAtomIdx()
                bond_a=mol.GetAtomWithIdx(mol_begin).GetSymbol()
                bond_b=mol.GetAtomWithIdx(mol_end).GetSymbol()
                #print(bond_a)
                #print(bond_b)
                if bond_a.lower() == 'c' and bond_b.lower() == 'o':
                    #print('find a desired one with C=O bond!')
                    continue
                elif bond_a.lower() == 'o' and bond_b.lower() == 'c':
                    #print('find a desired one with O=C bond!')
                    continue
                else:
                    #print("find a undesired one with", bond_a, bond_b)
                    return False
            else:
                #print('find a desired one!')
                continue
        return True
    

    # TODO: QED calculation
    def qed(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        qed = QED.qed(mol)
        return qed


    # TODO: Fsp3 calculation
    def fsp3(self, smiles):
        """
        Calculate the Fsp³ (fraction of sp³-hybridized carbons) for a given molecule.
        
        Parameters:
        - smiles (str): SMILES representation of the molecule.
        
        Returns:
        - fsp3 (float): Fraction of sp³-hybridized carbons in the molecule.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        
        # Calculate Fsp3 using RDKit's built-in function
        return rdMolDescriptors.CalcFractionCSP3(mol)


    # TODO: Fraction C-atoms calculation
    def fraction_c(self, smiles):
        """
        Calculate the fraction of carbon atoms in a molecule.
        
        Parameters:
        - smiles (str): SMILES representation of the molecule.
        
        Returns:
        - fraction_c (float): Fraction of carbon atoms in the molecule.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        
        total_atoms = mol.GetNumAtoms()
        if total_atoms == 0:
            return 0.0
        
        # Count the number of carbon atoms
        carbon_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        
        # Calculate fraction of carbon atoms
        return carbon_atoms / total_atoms
    

    # TODO: less than 1 7/8 membered ring, no ring greater than 8 membered
    def ring_size_check(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        one_7_8_membered_ring = 0
        zero_7_8_membered_ring = 0
        more_than_one_7_8_membered_ring = 0
        pass_check = True
        try:
                # Get the ring information
                ring_info = mol.GetRingInfo()

                # Get the ring sizes
                ring_sizes = []
                for ring in ring_info.AtomRings():
                    ring_sizes.append(len(ring))

                ring_sizes_set = sorted(ring_sizes)

                count_7_or_8 = sum(1 for element in ring_sizes_set if 7 <= element <= 8)
                count_greater_than_8 = sum(1 for element in ring_sizes_set if element > 8)

                if count_greater_than_8 == 0 and count_7_or_8 == 1:
                    one_7_8_membered_ring = 1

                elif count_greater_than_8 == 0 and count_7_or_8 < 1:
                    zero_7_8_membered_ring = 1

                else:
                    more_than_one_7_8_membered_ring = 1
                    pass_check = False

        except:
            more_than_one_7_8_membered_ring = 1
            pass_check = False

        return pass_check
    

    # TODO: Check undesired functional group
    def undesired_FG_check(self, smiles):
        pass_check = True
        
        filter1 = self.non_aromatic_double_bond_filter(smiles)
        filter2 = self.alogp(smiles)
        filter3 = self.if_NO_NN_in_non_aromatic_ring_or_acyclic(smiles)
        filter4 = self.if_contain_OCO(smiles)
        filter5 = self.if_contain_N3ring(smiles)
        filter6 = self.threeringcheck(smiles)
        filter7 = self.has_small_rings(smiles)
        filter8 = self.ring_size_check(smiles)
        filter9 = self.has_atom_in_three_rings(smiles)
        filter10 = self.divalent_nodes_fraction(smiles)

        if filter7 == True:
            pass_check = False
        elif filter10 <= 0.4:
            pass_check = False
        elif filter6 == False:
            pass_check = False
        elif filter8 == False:
            pass_check = False
        elif filter9 == True:
            pass_check = False
        elif filter1 == False:
            pass_check = False
        elif filter3 == True:
            pass_check = False
        elif filter4 == True:
            pass_check = False
        elif filter5 == True:
            pass_check = False
        elif filter2 < 0:
            pass_check = False

        return [pass_check, [filter1, filter2, filter3, filter4, filter5, filter6, filter7, filter8, filter9, filter10]]



    # TODO: Check undesired functional group
    def undesired_FG_details(self, FILE_PATH_READ, FILE_PATH_SAVE):
        df = pd.read_csv(FILE_PATH_READ, names=["SMILES", "Log Prob"], sep="\t")

        df['Filter1-10'] = df['SMILES'].apply(self.undesired_FG_check)

        df_failed = df[
            df['Filter1-10'].apply(
                lambda x: isinstance(x, list) and len(x) > 0 and x[0] == False
            )
        ].reset_index(drop=True)

        df_failed.to_csv(FILE_PATH_SAVE, sep='\t', header=False, index=False)
        print(f"File saved successfully to {FILE_PATH_SAVE}")

        return len(df_failed)