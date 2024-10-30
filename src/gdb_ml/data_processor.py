import pandas as pd
import ast
import re

COLUMN_NAME_INPUT = ["SMILES"]
COLUMN_NAME_OUTPUT = ["Predicted SMILES"]
COLUMN_NAME_TREAT = ["SMILES"]
SEPRATOR = '\t'
FILE_PATH_READ = "file_path_read"
FILE_PATH_SAVE = "file_path_save"


class DataProcessor:


    def __init__(self):
        # Initialize any instance variables here, if needed
        pass


    # TODO: Load the data
    def load_data(self, FILE_PATH_READ, SEPRATOR, file_type="csv"):
        """Loads data from the specified file type"""
        if file_type == "csv":
            # Read the file twice: once with header and once without
            df_with_header = pd.read_csv(FILE_PATH_READ, sep=SEPRATOR, nrows=1)
            df_without_header = pd.read_csv(FILE_PATH_READ, sep=SEPRATOR, header=None, nrows=1)

            # Check if column names appear as data in the first row without headers
            has_header = not df_with_header.columns.equals(df_without_header.iloc[0])

            if has_header:
                print("The file has a header row.")
                new_data = pd.read_csv(FILE_PATH_READ, sep=SEPRATOR)

            else:
                print("The file does not have a header row.")
                new_data = pd.read_csv(FILE_PATH_READ, names=COLUMN_NAME_INPUT, sep=SEPRATOR)
                
        elif file_type == "json":
            new_data = pd.read_json(FILE_PATH_READ)
        elif file_type == "excel":
            new_data = pd.read_excel(FILE_PATH_READ)
        else:
            raise ValueError("Unsupported file type. Supported types: csv, json, excel")
        return new_data


    # TODO: Load the file with correct number of lines when it may contain bad lines
    def load_file_with_badlines(self, FILE_PATH_READ):
        valid_lines = []
        total_lines = 0
        bad_lines = 0
    
        # Read the file using Python's open function
        with open(FILE_PATH_READ, 'r') as f:
            lines = f.readlines()
        
        print(f"Total rows in file: {len(lines)}")

        # Process each line (split by tab)
        data = [line.strip().split('\t') for line in lines]
        # Convert the processed data into a pandas DataFrame
        df = pd.DataFrame(data, columns = COLUMN_NAME_OUTPUT)
        return df
        # "['C', 'C', 'C', '1', 'C', '2', 'C', '3', 'C', 'C', '3', 'C', '1', 'C', '2', '(', 'C', '#', 'N', ')', 'C', '#',--> bad line"
        # df = pd.read_csv(FILE_PATH_READ, sep='\t', skip_blank_lines = False, error_bad_lines=True, warn_bad_lines=True)
        # print(df.shape) 


    # TODO: Remove the lines containing dots:
    def remove_dot_lines(df, COLUMN_NAME_TREAT):
        """
        Remove rows from the DataFrame where the specified column contains dots ('.').
        """
        # Find indices where the column contains a dot
        indices_to_remove = df[df[COLUMN_NAME_TREAT].str.contains(r'\.', na=False)].index
        # Drop those rows and return the resulting DataFrame
        df_cleaned = df.drop(indices_to_remove)
        print(f"Removed {len(indices_to_remove)} rows containing dots from '{COLUMN_NAME_TREAT}' column.")
        return df_cleaned


    # TODO: Remove the duplicates
    def remove_duplicates(df, COLUMN_NAME_TREAT, keep='last'):
        """
        Removes duplicate rows in the DataFrame.

        Parameters:
        keep (str, optional): Determines which duplicates (if any) to keep.
                            'first' keeps the first occurrence, 'last' keeps the last.
                            Default is 'last'.
        Returns:
        pd.DataFrame: A new DataFrame with duplicates removed.
        """
        # Drop duplicates based on the subset and keep specified rows
        df_cleaned = df.drop_duplicates(subset=COLUMN_NAME_TREAT, keep=keep).reset_index(drop=True)
        print(f"Removed duplicates based on columns: {COLUMN_NAME_TREAT}. Remaining rows: {len(df_cleaned)}")
        return df_cleaned


    # TODO: Tokenize the data
    def tokenize_smiles(self, smiles):
        # Define regular expression for SMILES tokens, atomwise
        pattern = r"(\[[^\]]*\]|Br|Cl|@@?|[=#$:\-\+\(\)/\\\.]|[A-Za-z0-9])"
        return re.findall(pattern, smiles)


    # TODO: Detokenize the data
    def detokenize_smiles(self, tokens_string):
        tokens_string_cleaned = tokens_string.strip('"')
        tokens_list = ast.literal_eval(tokens_string_cleaned)
        detokenized_string = ''.join(tokens_list)
        return detokenized_string


    # TODO: Save the output file
    def save_to_file(self, df, FILE_PATH_SAVE):
        """
        Saves the DataFrame to a specified file.

        Parameters:
        df (pd.DataFrame): The DataFrame to save.
        file_path (str): Path where the file will be saved.
        sep (str, optional): Separator for the output file. Default is '\t' (tab-separated).
        header (bool, optional): If False, no header is written. Default is False.
        index (bool, optional): If False, row indices are not written. Default is False.
        
        Returns:
        None
        """
        try:
            df.to_csv(FILE_PATH_SAVE, sep='\t', header=False, index=False)
            print(f"File saved successfully to {FILE_PATH_SAVE}")
        except Exception as e:
            print(f"An error occurred while saving the file: {e}")
