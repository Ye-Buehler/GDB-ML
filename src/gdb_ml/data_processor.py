import pandas as pd
import ast
import re

COLUMN_NAME = ["SMILES"]
SEPRATOR = '\t'
FILE_PATH_READ = "file_path_read"
FILE_PATH_READ = "file_path_save"

class DataProcessor:
    def __init__(self):
        # Initialize any instance variables here, if needed
        pass
    
    # TODO: Load the data
    def load_data(self, file_path, file_type="csv"):
        """Loads data from the specified file type"""
        if file_type == "csv":
            new_data = pd.read_csv(FILE_PATH_READ, names=COLUMN_NAME, sep=SEPRATOR)
        elif file_type == "json":
            new_data = pd.read_json(FILE_PATH_READ)
        elif file_type == "excel":
            new_data = pd.read_excel(FILE_PATH_READ)
        else:
            raise ValueError("Unsupported file type. Supported types: csv, json, excel")

    # TODO: Show the number of lines in the file, with or without the bad lines
    # "['C', 'C', 'C', '1', 'C', '2', 'C', '3', 'C', 'C', '3', 'C', '1', 'C', '2', '(', 'C', '#', 'N', ')', 'C', '#',--> bad line"""

    with open('/mnt/sdb/PROJECT_ML/ChEMBl_GDB13s_graph_learning/---TRANSFORMER_UBELIX---/test14-pure_ChEMBL_REDO_no_data_leakage/translation-beam50/pred_pure_ChEMBLa_51419.txt', 'r') as f:
        lines = f.readlines()
        print(f"Total rows: {len(lines)}")

    """df = pd.read_csv('TRANSFORMER_UBELIX/test5/translation-beam50/pred_pure_ChEMBL_53989.txt', sep='\t', skip_blank_lines = False, error_bad_lines=True, warn_bad_lines=True)
    print(df.shape) """


    # TODO: Remove the lines containing dots:
    """indices_to_remove = df1[df1['Source SMILES'].str.contains('\.', na=False)].index"""


    # TODO: Randomly sample the data
    #def random_sample(self, subset=None):


    # TODO: Remove the duplicates
    """def remove_duplicates(self, subset=None):
        """
        Removes duplicate rows in the data.
        :param subset: List of column names to check for duplicates.
                    If None, considers all columns.
        """
        self.data.drop_duplicates(subset=subset, inplace=True)

        df1 = df_output.drop_duplicates(subset=['Canonical SMILES'], keep='last').reset_index(drop=True)
        df1"""

    # TODO: Tokenize the data
    def tokenize_smiles(smiles):
        # Define regular expression for SMILES tokens, atomwise
        pattern = r"(\[[^\]]*\]|Br|Cl|@@?|[=#$:\-\+\(\)/\\\.]|[A-Za-z0-9])"
        return re.findall(pattern, smiles)

    # TODO: Detokenize the data
    def detokenize_smiles(tokens_string):
        tokens_string_cleaned = tokens_string.strip('"')
        tokens_list = ast.literal_eval(tokens_string_cleaned)
        detokenized_string = ''.join(tokens_list)
        return detokenized_string

    # TODO: Save the output file

    """    df_ligands.to_csv(
        FILE_PATH_READ",
            sep='\t', header = False, index=False)"""