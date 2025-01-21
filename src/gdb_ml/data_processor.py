import pandas as pd
from rdkit import Chem
import ast
import re
import sys
import os
from tqdm import tqdm
from pandarallel import pandarallel

# Initialize pandarallel
pandarallel.initialize(progress_bar=True)  # Enable progress bar for better tracking

FILE_PATH_READ = "file_path_read"
FILE_PATH_SAVE = "file_path_save"


class DataProcessor:


    def __init__(self):
        # Initialize any instance variables here, if needed
        pass


    # TODO: Load the data
    def load_data(self, FILE_PATH_READ, SEPRATOR='\t', file_type="csv", has_header=0, add_header=1, COLUMN_NAME_INPUT= ["SMILES"]):
        """Loads data from the specified file type"""
        if file_type == "csv":
            # Read the file twice: once with header and once without
            df_with_header = pd.read_csv(FILE_PATH_READ, sep=SEPRATOR, nrows=1)
            df_without_header = pd.read_csv(FILE_PATH_READ, sep=SEPRATOR, header=None, nrows=1)

            if has_header==1:
                print("The file has a header row.")
                new_data = pd.read_csv(FILE_PATH_READ, sep=SEPRATOR)

            else:
                print("The file does not have a header row.")
                if add_header==0:
                    new_data = pd.read_csv(FILE_PATH_READ, sep=SEPRATOR)
                    print("No header has been added.")
                else:
                    new_data = pd.read_csv(FILE_PATH_READ, names=COLUMN_NAME_INPUT, sep=SEPRATOR)
                    print("A header has been added.")
                
        elif file_type == "json":
            new_data = pd.read_json(FILE_PATH_READ)
        elif file_type == "excel":
            new_data = pd.read_excel(FILE_PATH_READ)
        elif file_type == "sdf":
            # Load the molecules from the SDF file
            suppl = Chem.SDMolSupplier(FILE_PATH_READ)
            # Extract properties into a DataFrame
            data = []
            for mol in suppl:
                if mol is not None:
                    props = mol.GetPropsAsDict()
                    data.append(props)
            new_data = pd.DataFrame(data)
        else:
            raise ValueError("Unsupported file type. Supported types: csv, json, excel, sdf")
        return new_data


    # TODO: Load the file with correct number of lines when it may contain bad lines
    def load_file_with_badlines(self, FILE_PATH_READ, COLUMN_NAME_OUTPUT):
    
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
        Saves the DataFrame to a csv file.
        """
        try:
            df.to_csv(FILE_PATH_SAVE, sep='\t', header=False, index=False)
            print(f"File saved successfully to {FILE_PATH_SAVE}")
        except Exception as e:
            print(f"An error occurred while saving the file: {e}")


    # TODO: Save a list into a txt file
    def save_list(self, list, FILE_PATH_SAVE):
        # Write the list to a text file, each element on a new line
        with open(FILE_PATH_SAVE, 'w') as file:
            for item in list:
                file.write(f"{item}\n")
        print(f"File saved successfully to {FILE_PATH_SAVE}")


    # TODO: Read a list from a txt file
    def read_list(self, FILE_PATH_READ):
        # Reading the list from a text file (line by line)
        with open(FILE_PATH_READ, 'r') as file:
            list = [line.strip() for line in file]  # Remove any leading/trailing whitespace
        return list
    

    # TODO: append all the dataframes from a folder
    def append_dfs_in_folder(self, FOLDER_PATH, COLUMN_NAME_INPUT, SEPRATOR='\t', FILETYPE = ".txt") -> pd.DataFrame:

        iteration = 0
        row_count = 0

        for filename in os.listdir(FOLDER_PATH):
            if filename.endswith(FILETYPE):
                
                iteration += 1
                if iteration == 1:
                    df = pd.read_csv(FOLDER_PATH + filename, sep=SEPRATOR, names=COLUMN_NAME_INPUT)
                    #print("basis lenght = " + str(len(df)))
                    row_count += len(df)

                if iteration != 1:
                    df2 = pd.read_csv(FOLDER_PATH + filename, sep=SEPRATOR,names=COLUMN_NAME_INPUT)
                    #print("added lenght = " + str(len(df2)))
                    df = pd.concat([df, df2], axis=0, ignore_index=True)
                    row_count += len(df2)

        print("Total enteries that have been merged: \t" + str(row_count))
        return df
    

    # TODO: detokenize the files and append the log prob
    def detokenize_append_log_prob(self, FILE_PATH_READ_1, FILE_PATH_READ_2, FILE_PATH_SAVE):
        # TODO: Compute the metric validity and canonicalize the SMILES
        def validity(df):

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
        
        COLUMN_NAME_OUTPUT = ["SMILES"]
        df = self.load_file_with_badlines(FILE_PATH_READ_1, COLUMN_NAME_OUTPUT)
        df['Detokenized']=""
        count = 0
        my_list_invaild = []
        for q in tqdm(range(0, len(df)), desc = 'Loop 1'):
            try:  
                tokens_string = df['SMILES'][q]
                detokenized_string = self.detokenize_smiles(tokens_string)
                df.loc[q,'Detokenized'] = detokenized_string 
            except:
                count += 1
                my_list_invaild.append(q)
                continue
        df_cs = self.load_data(FILE_PATH_READ_2, COLUMN_NAME_INPUT=['log prob'])
        df = pd.DataFrame({'SMILES': df['Detokenized'], 'Log Probs': df_cs['log prob']})
        validity, df_valid = validity(df)
        self.save_to_file(df_valid, FILE_PATH_SAVE)
        return validity, df_valid
    
    # TODO: combine the log prob values and count correspondingly
    def log_prob_count(self, FOLDER_PATH):

        # Function to process each file
        def process_file(file_path):
            # Read the file into a DataFrame (assuming tab-separated)
            df = pd.read_csv(file_path, sep='\t', header=None, usecols=[1], names=["Log Prob"])

            # Round the "Log Prob" to 1 decimal place
            df['Log Prob'] = df['Log Prob'].round(1)

            # Count occurrences of each log prob value
            return df['Log Prob'].value_counts().reset_index(name='Count').rename(columns={'index': 'Log Prob Category'})

        # Get a list of all txt files in the folder
        files = [os.path.join(FOLDER_PATH, f) for f in os.listdir(FOLDER_PATH) if f.endswith('.txt')]

        # Process files in parallel using pandarallel and convert the result to a list
        results = pd.Series(files).parallel_map(process_file).tolist()

        # Concatenate the results into a single DataFrame
        combined_results = pd.concat(results)

        # Group by "Log Prob Category" and sum the counts
        final_result = combined_results.groupby('Log Prob Category', as_index=False)['Count'].sum()

        # Sort by "Log Prob Category"
        final_result = final_result.sort_values(by='Log Prob Category')

        return final_result