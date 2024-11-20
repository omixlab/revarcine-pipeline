import os
import sys
import pandas as pd
import configparser
import warnings as wn

wn.filterwarnings('ignore')
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

from exe.models import *
from exe.pre_process import pre_processing
from exe.pos_process import pos_processing

def revarcine(PARAMS):


    print("                                                            ")
    print("    ____  ____ __     ___    ____   ____ ___ _   _ _____    ")
    print("   |  _ \| ____\ \   / / \  |  _ \/ ___|_ _| \ | | ____|   ")
    print("   | |_) |  _|  \ \ / / _ \ | |_) | |    | ||  \| |  _|     ")
    print("   |  _ <| |___  \ V / ___ \|  _ <| |___ | || |\  | |___    ")
    print("   |_| \_\_____|  \_/_/   \_\_| \_\\____|___|_| \_|_____|   ")
    print("                                                            ")
    print("____________________________________________________________")
    print("           Revarcine - Reverse Vaccinology Tool v1          ")
    print("                                                            ")

    # Set Variables
    INPUT_NAME  = PARAMS[0] # .fasta file
    INPUT_GRAM  = PARAMS[1] # G+ ou G-python -W ignore revarcine.py data.fasta 'G+' '110'
    ARGPRED     = PARAMS[2] # 100 -- Analysis combination --> Default: All Analysis -> check the documentation

    models_gp   = ['spp','cgp', 'scp']
    models_gn   = ['spgn','cgn', 'scn', 'bbn', 'ahn']
    list_inter  = list()
    analysis    = list()
    results     = list()
    argms       = ARGPRED
    select_gram = INPUT_GRAM
    LIMITE      = 5
    

    if argms in ('000', '00000'):
        print('Warning !')
        print(f'You did not selected any analysis, the scirpt is in Default Mode: {argms}')
        sys.exit()
    else:
        pass
    
    if select_gram == 'G+':
        models_list = models_gp

    else:
        models_list = models_gn

    # Import the input (.fasta file) and Run pre-processing
    try:
        path  = 'data/input/' + INPUT_NAME
        fasta = open(path)
        try:
            df_pre = pre_processing(fasta)

        except Exception as exp:
            print(exp)
            print("Unnable to pre-process the data")
            sys.exit()

    except Exception as exc:
        print(exc)
        print("Unnable to load the .fasta file")
        sys.exit()

    # Split DataFrame into smaller chunks
    def split_dataframe(df, chunk_size):
        return [df.iloc[i:i+chunk_size] for i in range(0, df.shape[0], chunk_size)]
    
    # Revarcine G+
    class Revar_Pos:
        def cgp(df):
            try:
                sgp_result = models.cgp(df)            

            except Exception as exp:
                print(exp)
                print("Unable to execute Signal Peptide Cut Position model")
                sys.exit()
            return sgp_result
        
        def scp(df):
            try:
                scp_result = models.scp(df)
                
            except Exception as exp:
                print(exp)
                print("unable to execute Subcellular Location model")
                sys.exit()
            return scp_result
        
        def spp(df):
            try:
               spp_result = models.spp(df)

            except Exception as exp:
                print(exp)
                print("Unable to execute Signal Peptide model")
                sys.exit()
            return spp_result
        
    # Revarcine G-
    class Revar_Neg:

        def cgn(df):
            try:
               cgn_result =  models.cgn(df)
                
            except Exception as exp:
                print(exp)
                print("unable to execute Signal Peptide Cut Position model")
                sys.exit()
            return cgn_result
        
        def scn(df):
            try:
               scn_result = models.scn(df)

            except Exception as exp:
                print(exp)
                print("Unable to execute Subcellular Location model")
                sys.exit()
            return scn_result
        
        def spgn(df):
            try:
                spn_result = models.spgn(df)
                
            except Exception as exp:
                print(exp)
                print("unable to execute Signal Peptide model")
                sys.exit()
            return spn_result

        def bbn(df):
            try:
                bbn_result = models.bbn(df)
                
            except Exception as exp:
                print(exp)
                print("unable to execute Beta Berral model")
                sys.exit()
            return bbn_result
            
        def ahn(df):
            try:
                ahn_result = models.ahn(df)
                
            except Exception as exp:
                print(exp)
                print("unable to execute Alpha Helix model")
                sys.exit()
            return ahn_result
        
    
    # Split the args in a list
    list_model = [f for f in argms]

    # For next update
    if len(list_model)>LIMITE:
        print('\nWARNING:')
        print(f'The knew Model will be integrated sooner, please select another configuration than {list_model}')
        sys.exit()
    else:
        # Get the index of models required
        for i in list_model:
            if i == "1":
                element_index             = list_model.index(i)
                # Replace the value of the index by zero to avoid "infinty loop or false positive"
                list_model[element_index] = "0" 
                list_inter.append(element_index)

        # Get the names of required models
        for j in list_inter:
            analysis.append(models_list[j])

        
        # Get a specific method in the class and execute
        print('\nRunning the models')
        print('_____________________')
        
        if select_gram == "G+" and len(list_model)<LIMITE:
            # Going through each sequence and execute the models
            for k in analysis:
                
                print(f'\nThe model {k} is processing') 
                out_model = getattr(Revar_Pos, k)
                # Check point for the result
                if callable(out_model):
                    results.append(out_model(df_pre))
                else:
                    print(f'{out_model} is not callable or does not exist')
                
        elif select_gram == "G-"and len(list_model)==LIMITE:
            # Processing each smaller DataFrame
            for k in analysis:
                print('\n')
                print(f'The model {k} is processing')              
                out_model = getattr(Revar_Neg, k)
                if callable(out_model):
                    results.append(out_model(df_pre))
                else:
                    print(f'{out_model} is not callable or does not exist')
            #print(results)
        else:
            print("Please provide the right gram classification")
            print("Remember that for +G, the sequency must have only three positions and five for G-. Eg: 101,10001")
            sys.exit()
    
    # Pos-processing the result
    try:    
        pos_processing(df_pre ,results, analysis)
    except Exception as exp:
        print(exp)
        print("Unnable to pos-process the model result")
        sys.exit()
  

# Running the main function
if str(sys.argv[0]) == "revarcine.py":
    PARAMS = sys.argv[1:]
    
    revarcine(PARAMS)

print('DONE!')
