import pandas as pd
import numpy as np
import json
import os

MOLS_PER_GRAPH = 10
TOTAL_DATAPOINTS = 2000000

class GraphMapping:


    def __init__(self):
        # Initialize any instance variables here, if needed
        pass

    # TODO: Split the datapoints into three sets
    def datapoints_split(self, df, TOTAL_DATAPOINTS, MOLS_PER_GRAPH):

        df["Number of Datapoints"]= MOLS_PER_GRAPH
        df['Accumulated Sum of Datapoints'] = df['Number of Datapoints'].cumsum().astype(int)
        df = df[0:int(TOTAL_DATAPOINTS/MOLS_PER_GRAPH)]

        # Create a new column 'Set Assigned' with repeated values 'train', 'val', 'test'
        sets = ['train','train','train','train','train','train','train','train','train','train','train','train','train','train','train','train','train', 'val', 'val','test']
        df['Set Assigned'] = [sets[i % len(sets)] for i in range(len(df))]

        # Split the DataFrame into three based on the 'Set Assigned' column
        df_train = df[df['Set Assigned'] == 'train'].reset_index(drop=True)
        df_train = df_train.sort_values(by='Number of Values', ascending=False).reset_index(drop = True)
        df_train['Accumulated Sum of Datapoints'] = df_train['Number of Datapoints'].cumsum().astype(int)
        df_train['Accumulated Sum'] = df_train['Number of Values'].cumsum().astype(int)
        
        df_val = df[df['Set Assigned'] == 'val'].reset_index(drop=True)
        df_val['Accumulated Sum'] = df_val['Number of Values'].cumsum().astype(int)
        df_val = df_val.sort_values(by='Number of Values', ascending=False).reset_index(drop = True)
        df_val['Accumulated Sum of Datapoints'] = df_val['Number of Datapoints'].cumsum().astype(int)
        df_val['Accumulated Sum'] = df_val['Number of Values'].cumsum().astype(int)

        df_test = df[df['Set Assigned'] == 'test'].reset_index(drop=True)
        df_test['Accumulated Sum'] = df_test['Number of Values'].cumsum().astype(int)
        df_test = df_test.sort_values(by='Number of Values', ascending=False).reset_index(drop = True)
        df_test['Accumulated Sum of Datapoints'] = df_test['Number of Datapoints'].cumsum().astype(int)
        df_test['Accumulated Sum'] = df_test['Number of Values'].cumsum().astype(int)

        df_train['Actual Number of Datapoints'] = ""
        df_train['Actual cumulated Datapoints'] = ""
        df_train['Actual Number of Datapoints'] = np.where(df_train['Number of Values'] > MOLS_PER_GRAPH, MOLS_PER_GRAPH, df_train['Number of Values'])
        df_train['Actual cumulated Datapoints'] = df_train['Actual Number of Datapoints'].cumsum().astype(int)

        df_val['Actual Number of Datapoints'] = ""
        df_val['Actual cumulated Datapoints'] = ""
        df_val['Actual Number of Datapoints'] = np.where(df_val['Number of Values'] > MOLS_PER_GRAPH, MOLS_PER_GRAPH, df_val['Number of Values'])
        df_val['Actual cumulated Datapoints'] = df_val['Actual Number of Datapoints'].cumsum().astype(int)

        df_test['Actual Number of Datapoints'] = ""
        df_test['Actual cumulated Datapoints'] = ""
        df_test['Actual Number of Datapoints'] = np.where(df_test['Number of Values'] > MOLS_PER_GRAPH, MOLS_PER_GRAPH, df_test['Number of Values'])
        df_test['Actual cumulated Datapoints'] = df_test['Actual Number of Datapoints'].cumsum().astype(int)
        
        list_train = df_train['Key'].tolist()
        list_val = df_val['Key'].tolist()
        list_test = df_test['Key'].tolist()

        return [df_train, df_val, df_test, list_train, list_val, list_test]

    
    # TODO: Split the datapoints into three sets
    def check_mols_from_graph(self, keys_to_merge, FOLDER_PATH, MOLS_PER_GRAPH):

        merged_dict = {key: [] for key in keys_to_merge}
        dicts = []

        for filename in os.listdir(FOLDER_PATH):
            with open(FOLDER_PATH + filename, 'r') as file:
                loaded_dict = json.load(file)
                dicts.append(loaded_dict)

        for d in dicts:
            for key, value in d.items():
                if key in merged_dict:
                    # If key already exists, append the new value to the list
                    merged_dict[key].append(value)
                else:
                    # If key doesn't exist, add it with its value
                    merged_dict[key] = value       

        # Flatten the lists for each key in merged_values
        flattened_values = {key: [item for sublist in merged_dict[key] for item in sublist] for key in merged_dict}
        flattened_values

        # Truncate the values in the dictionary
        truncated_dict = {key: values[:MOLS_PER_GRAPH] for key, values in flattened_values.items()}
        truncated_dict

        # Initialize variables to keep track of the accumulated sum
        accumulated_sum = 0

        # Get the number of unique keys
        num_unique_keys = len(truncated_dict)

        # Iterate through the dictionary to count values and compute the accumulated sum
        for key, value in truncated_dict.items():
            num_values = len(value)  # Number of values for the current key
            accumulated_sum += num_values  # Add to the accumulated sum
            print(f"Key: {key}, Number of Values: {num_values}, Accumulated Sum: {accumulated_sum}")     
        print(f"Number of Unique Keys: {num_unique_keys}")
        return truncated_dict