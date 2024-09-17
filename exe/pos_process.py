import os
import numpy as np
import pandas as pd
import warnings as wn

wn.filterwarnings('ignore')
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

def pos_processing(df_pre ,results ,analysis):
    print("\nStart the Pos-Processing")
    
    dt             = df_pre.copy()
    names          = list() 
    
    model_names    = {'spp':'Positive Signal Pepitede',
                    'cgp':'Gram Positive Cut',
                    'scp':'Positive Suncellular Location',
                    'spgn':'Negative Signal Pepitede',
                    'cgn':'Gram Negative Cut',
                    'scn':'Negative Suncellular Location',
                    'bbn':'Beta Barrel'}
    
    scn_categories = {
        0:'Cell inner membrane',
        1:'Cell outer membrane',
        2:'Cytoplasm',
        3:'Periplasm',
        4:'Secreted',
        5: 'Other'
    }
    
    scp_categories = {
        0:'Cell inner membrane',
        1:'Cell outer membrane',
        2:'Cytoplasm',
        3:'Membrane',
        4:'Secreted',
        5:'Periplasm',
        6:'Other'
    }
    
    # Get the model full names
    for symbol in analysis:
        names.append(model_names[symbol])
    
    # Processar para remover listas internas com um único elemento
    processed_data = []
    for sublist in results:
        processed_data.append([item[0] if isinstance(item, list) and len(item) == 1 else item for item in sublist])

    # Criar DataFrame e transpor
    df             = pd.DataFrame(processed_data).transpose()
    df.columns = analysis
    
    
    # check if the code runned spp, spn,scp or scn model
    # For the ajustment and classification
    
    if any(col in analysis for col in ["spp", "spgn", "scp", "scn"]):
        if "spp" in analysis and "cgp" in analysis:
            df.loc[df['spp']  == 0, 'cgp'] = np.nan
        if "spgn" in analysis and "cgn" in analysis:
            df.loc[df['spgn'] == 0, 'cgn'] = np.nan
        if "scp" in analysis: 
            # Substituindo os valores numéricos pelas referências categóricas
            df['scp'] = df['scp'].replace(scp_categories)
        if "scn" in analysis:
            # Substituindo os valores numéricos pelas referências categóricas
            df['scn'] = df['scn'].replace(scn_categories)
    
    # Rename the columns
    df.columns = names
    
    # Reset the index of df_results to make the index a column for merging
    df = df.reset_index().rename(columns={'index': 'reference'})
    
    dt['reference'] = dt['reference'].apply(lambda x: int(x.split('_')[-1]) - 1)

    # Merge df1 and df2 based on the "reference" column
    df_result = pd.merge(dt, df, on='reference').drop(columns=['reference'])
    
    # Removendo os espaços da coluna "sequencia"
    df_result['sequence'] = df_result['sequence'].str.replace(' ', '')

    
    print(df_result.head())

    # Verifica se o diretório existe, caso contrário, cria-o
    # antes de salvar o resultado
    output_dir = 'data/output'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    df_result.to_csv(output_dir + '/Results.csv',index = False)
    
    print("Pos-Processing is done")
    
if __name__ == "__main__":
    pos_processing()
