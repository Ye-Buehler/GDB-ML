import pandas as pd
import numpy as np


MOLS_PER_GRAPH = 10

class GraphMapping:


    def __init__(self):
        # Initialize any instance variables here, if needed
        pass


    def datapoints_split(self, df, MOLS_PER_GRAPH):

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
        # Add new column 'B' where values from 'A' are capped at 50
        df_val['Actual Number of Datapoints'] = np.where(df_val['Number of Values'] > MOLS_PER_GRAPH, MOLS_PER_GRAPH, df_val['Number of Values'])
        df_val['Actual cumulated Datapoints'] = df_val['Actual Number of Datapoints'].cumsum().astype(int)

        df_test['Actual Number of Datapoints'] = ""
        df_test['Actual cumulated Datapoints'] = ""
        # Add new column 'B' where values from 'A' are capped at 50
        df_test['Actual Number of Datapoints'] = np.where(df_test['Number of Values'] > MOLS_PER_GRAPH, MOLS_PER_GRAPH, df_test['Number of Values'])
        df_test['Actual cumulated Datapoints'] = df_test['Actual Number of Datapoints'].cumsum().astype(int)
        
        list_train = df_train['Key'].tolist()
        list_val = df_val['Key'].tolist()
        list_test = df_test['Key'].tolist()

        return {
            "train dataframe": df_train,
            "val dataframe ": df_val,
            "test dataframe": df_test,
            "train graph list": list_train,
            "val graph list": list_val,
            "test graph list": list_test
            }