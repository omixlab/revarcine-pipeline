import pickle
import numpy as np
import pandas as pd
import configparser
import seaborn as sns
import tensorflow as tf
import matplotlib.pyplot as plt


from tensorflow import keras
from tensorflow.keras import layers

from sklearn import metrics
from sklearn.utils import shuffle
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, matthews_corrcoef, f1_score

from keras.callbacks import ModelCheckpoint
from keras.preprocessing.text import Tokenizer
from keras.models import Sequential, save_model
from keras_preprocessing.sequence import pad_sequences
from keras.layers import Dense, Embedding, GlobalMaxPooling1D, Dropout, LSTM


class TransformerBlock(layers.Layer):
    def __init__(self, embed_dim, num_heads, ff_dim, rate=0.1):
        super(TransformerBlock, self).__init__()
        self.embed_dim = embed_dim
        self.num_heads = num_heads
        self.ff_dim = ff_dim
        self.rate = rate
        self.att = layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim)
        self.ffn = keras.Sequential(
            [layers.Dense(ff_dim, activation="relu"), layers.Dense(embed_dim),]
        )
        self.layernorm1 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = layers.LayerNormalization(epsilon=1e-6)
        self.dropout1 = layers.Dropout(rate)
        self.dropout2 = layers.Dropout(rate)

    def call(self, inputs, training):
        attn_output = self.att(inputs, inputs)
        attn_output = self.dropout1(attn_output, training=training)
        out1 = self.layernorm1(inputs + attn_output)
        ffn_output = self.ffn(out1)
        ffn_output = self.dropout2(ffn_output, training=training)
        return self.layernorm2(out1 + ffn_output)
    
    def get_config(self):

        config = super().get_config().copy()
        config.update({
            'embed_dim': self.embed_dim,
            'num_heads': self.num_heads,
            'ff_dim': self.ff_dim,
            'rate': self.rate,
        })
        return config

def run():

     # Get information from config.txt file
    arq = configparser.RawConfigParser()
    arq.read('exe/config.txt')

    file_path = arq.get('FILE_SP_GP','file_path')

    # hyperparameters
    SAMPLES_PER_GROUP = int(arq.get('CONSTANTS_SP_GP','SAMPLES_PER_GROUP'))
    MAX_VOCAB_SIZE = int(arq.get('CONSTANTS_SP_GP','MAX_VOCAB_SIZE'))
    SEQUENCE_SIZE = int(arq.get('CONSTANTS_SP_GP','SEQUENCE_SIZE'))
    EMBEDDING_SIZE = int(arq.get('CONSTANTS_SP_GP','EMBEDDING_SIZE'))

    def data_process():

        # Data cleaning
        df = pd.read_csv(file_path, index_col=False)
        df = df.drop(['Unnamed: 0'], axis=1)

        classes_zero = df[df['signal_peptide'] == 0]
        classes_one  = df[df['signal_peptide'] == 1]

        final_dataframe = pd.concat(
        [
            df[df['signal_peptide']==0].sample(SAMPLES_PER_GROUP),
            df[df['signal_peptide']==1].sample(SAMPLES_PER_GROUP)
        ]
        )

        final_dataframe = shuffle(
            final_dataframe
        )

        # split the data
        df_train, df_test = train_test_split(final_dataframe, test_size=0.25, random_state=60)

        df_test['sequence'] = df_test['sequence'].str.replace(" ","")
        df_test.to_csv('gram_positive_SP.csv')

        return df_train, df_test
    
    def MODEL():

        df_train, df_test = data_process()

        # x train tokenizer
        tokenizer = Tokenizer(num_words=MAX_VOCAB_SIZE, char_level=True)
        tokenizer.fit_on_texts(df_train['sequence'])
        sequences = tokenizer.texts_to_sequences(df_train['sequence'])
        x_train = pad_sequences(sequences, maxlen=SEQUENCE_SIZE, padding='post', truncating='post', value=0)

        # x test tokenizer
        sequences = tokenizer.texts_to_sequences(df_test['sequence'])
        x_test = pad_sequences(sequences, maxlen=SEQUENCE_SIZE, padding='post', truncating='post', value=0)

        # y to numpy array
        Y_train = df_train['signal_peptide'].to_numpy()
        y_test = df_test['signal_peptide'].to_numpy()


        embed_dim = x_train.shape[1]  # Embedding size for each token
        num_heads = 1  # Number of attention heads
        ff_dim = 32  # Hidden layer size in feed forward network inside transformer

        model = Sequential()
        model.add(Embedding(input_dim=MAX_VOCAB_SIZE, output_dim=EMBEDDING_SIZE, input_length=SEQUENCE_SIZE, trainable = True))
        model.add(TransformerBlock(embed_dim, num_heads, ff_dim))
        model.add(Dropout(0.1))
        model.add(LSTM(50))
        model.add(Dense(100, input_dim=embed_dim, activation='relu'))
        model.add(Dropout(0.1))
        model.add(Dense(units = 1, activation = 'sigmoid'))
        model.summary()

        model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['acc'])

        model.fit(x_train, Y_train, validation_data = (x_test, y_test), batch_size = 41, epochs = 100)

        #model.load_weights('weights.best.hdf5')
        print(model.evaluate(x_test, y_test))

        # save tokenizer
        with open('spp.pkl', 'wb') as writer:
            writer.write(pickle.dumps(tokenizer))
        
        return model, x_test, y_test
    
    model, x_test, y_test = MODEL()
    model.save("spp.h5")

    # Prediction
    predict = model.predict(x_test)

    # filter
    #predict = predict >= 0.5

    #predict = pd.DataFrame(predict)
    #predict = predict.rename(columns={0: 'Peptide'})

    # Correction
    #matthews_corrcoef(y_test, predict)

    #acuracia = round(metrics.accuracy_score(y_test, predict),3)

    #print('\n--------------------------------------------------------\n')
    #print(predict.head(5))
    #print('\n')
    #print(f'Acurácia do model é de {acuracia}')
    #print(f"E o valor de f1_score é de {f1_score(y_test, predict, average='macro')}")
# ------------------------- Run the code -----------------------------------

if __name__ == '__main__':
    run()
