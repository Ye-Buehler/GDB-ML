from rdkit import Chem

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