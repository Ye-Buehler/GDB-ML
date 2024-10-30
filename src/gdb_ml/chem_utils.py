from rdkit import Chem
import re

class ChemUtils:
    def __init__(self):
        # Initialize any instance variables here, if needed
        pass

    # TODO: Canonicalize molecules
    def canonicalize_smiles(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string provided")
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        return canonical_smiles

    # TODO: Neutralizing molecules
    def neutralize_atoms(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
        at_matches = mol.GetSubstructMatches(pattern)
        at_matches_list = [y[0] for y in at_matches]
        if len(at_matches_list) > 0:
            for at_idx in at_matches_list:
                atom = mol.GetAtomWithIdx(at_idx)
                chg = atom.GetFormalCharge()
                hcount = atom.GetTotalNumHs()
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(hcount - chg)
                atom.UpdatePropertyCache()
        neutralized_smiles = Chem.MolToSmiles(mol)
        return neutralized_smiles

    # TODO: Convert the molecules into character-based graphs
    def graph_convert(self, smiles):
        updated_smiles = smiles.replace('.[Li+]', '')
        updated_smiles = updated_smiles.replace('.[Na+]', '')
        updated_smiles = updated_smiles.replace('.[K+]', '')
        updated_smiles = updated_smiles.replace('.[F-]', '')
        updated_smiles = updated_smiles.replace('.[Cl-]', '')
        updated_smiles = updated_smiles.replace('.[Br-]', '')
        updated_smiles = updated_smiles.replace('.[I-]', '')
        #updated_smiles = updated_smiles.replace('F', 'C')
        updated_smiles = updated_smiles.replace('Cl', 'C')
        updated_smiles = updated_smiles.replace('Br', 'C')
        #updated_smiles = updated_smiles.replace('I', 'C')
        #updated_smiles = updated_smiles.replace('N', 'C')
        #updated_smiles = updated_smiles.replace('O', 'C')
        #updated_smiles = updated_smiles.replace('P', 'C')
        #updated_smiles = updated_smiles.replace('S', 'C')
        #updated_smiles = updated_smiles.replace('c', 'C')
        #updated_smiles = updated_smiles.replace('n', 'C')
        #updated_smiles = updated_smiles.replace('o', 'C')
        updated_smiles = re.sub('\[.*?\]', 'C', updated_smiles)
        updated_smiles = re.sub('[A-Za-z]', 'C', updated_smiles)
        updated_smiles = updated_smiles.replace('-', '')
        updated_smiles = updated_smiles.replace('=', '')
        character_based_smiles = updated_smiles.replace('#', '')
        return character_based_smiles