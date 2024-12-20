from Bio import SeqIO
import re
import os
import subprocess
import pandas as pd
import warnings as wn

wn.filterwarnings('ignore')

def pre_processing(fasta):

    print("Start the Pre-Processing")
    
    # Baixar outros modelos do Google Drive: URL da pasta do Google Drive
    folder_url = 'https://drive.google.com/drive/u/0/folders/1g2pO1Kbtp1MRNSf2cLlVUnGhGSLPIts1'
    output_dir = 'exe/model'
    
    # Verifica se o diretório existe, caso contrário, cria-o
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Verifica se o diretório já contém arquivos .keras
    keras_files = [f for f in os.listdir(output_dir) if f.endswith('.keras')]

    if keras_files:
        print(f"The directory already contains .keras files. Downloading will not be necessary.\n")
    else:
        # Comando para baixar a pasta usando gdown
        command = ['gdown', '--folder', folder_url, '-O', output_dir]

        # Executa o comando
        try:
            subprocess.run(command, check=True)
            print(f"Download completed successfully. Files saved in: {output_dir}")
        except subprocess.CalledProcessError as e:
            print(f"Error executing command: {e}")  
    
    # Preparing the data.fasta

    records = list(SeqIO.parse(fasta, 'fasta'))
    
    df             = pd.DataFrame(columns=['id', 'sequence'])
    regex          = re.compile('>')
    sequence_id    = [record.id for record in records]
    sequence       = [str(record.seq) for record in records]
    # Convert all strings to lowercase
    sequence = [s.lower() for s in sequence]
    # Add a space between each character of the strings
    sequence = [' '.join(s) for s in sequence]

    df['id']       = sequence_id
    df['sequence'] = sequence
    # Adding the "reference" column
    df['reference'] = 'Sequence_' + (df.index + 1).astype(str)
    
    df.to_csv('data/input/df.csv', index = False) 
    
    print("Pre-Processing is done")

    return df

if __name__ == "__main__":
    pre_processing()
