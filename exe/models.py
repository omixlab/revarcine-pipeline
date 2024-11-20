import os
import pickle
import numpy as np
import configparser
import warnings as wn

wn.filterwarnings('ignore')
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

from tensorflow import keras
from keras.models import load_model
from tensorflow.keras import layers
from keras_preprocessing.sequence import pad_sequences


class TransformerBlock(layers.Layer):
    def __init__(self, embed_dim, num_heads, ff_dim, rate=0.1, **kwargs):
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


class models():

    def cgn(df): # gram negative cut position function model
        model = keras.models.load_model(
            "exe/model/cgn.keras", custom_objects={"TransformerBlock": TransformerBlock}
            ) # Loading model 

        with open('exe/tokenizer/cgn.pkl', 'rb') as reader:
            tokenizer = pickle.loads(reader.read()) # Loading tokenizer

        SEQUENCE_SIZE  = 50
        
        # tokenizer
        sequences = tokenizer.texts_to_sequences(df['sequence'])
        x         = pad_sequences(sequences, maxlen=SEQUENCE_SIZE, padding='post', truncating='post', value=0)
        # predict
        predict = model.predict(x) 
        result  = np.argmax(predict,axis=1).tolist()
        return result
    

    def cgp(df): #cut position function model
        
        
        model = keras.models.load_model(
            "exe/model/cgp.keras", custom_objects={"TransformerBlock": TransformerBlock}
            ) # Loading model 

        with open('exe/tokenizer/cgp.pkl', 'rb') as reader:
            tokenizer = pickle.loads(reader.read()) # Loading tokenizer

        SEQUENCE_SIZE  = 50

        # tokenizer
        sequences = tokenizer.texts_to_sequences(df['sequence'])
        x         = pad_sequences(sequences, maxlen=SEQUENCE_SIZE, padding='post', truncating='post', value=0)

        # predict
        predict = model.predict(x) 
        result  = np.argmax(predict,axis=1).tolist()
        return result

    def scn(df): # gram negative suncellular location function model
        
        model = keras.models.load_model(
            "exe/model/scn.keras", custom_objects={"TransformerBlock": TransformerBlock}
            ) # Loading model 

        with open('exe/tokenizer/scn.pkl', 'rb') as reader:
            tokenizer = pickle.loads(reader.read()) # Loading tokenizer
    

        SEQUENCE_SIZE  = 100
        
        # tokenizer
        sequences = tokenizer.texts_to_sequences(df['sequence'])
        x         = pad_sequences(sequences, maxlen=SEQUENCE_SIZE, padding='post', truncating='post', value=0)

        # predict
        predict = model.predict(x) 
        result  = np.argmax(predict,axis=1).tolist()
        return result

    def scp(df): # gram positive suncellular location function model
        
        model = keras.models.load_model(
            "exe/model/scp.keras", custom_objects={"TransformerBlock": TransformerBlock}
            ) # Loading model 

        with open('exe/tokenizer/scp.pkl', 'rb') as reader:
            tokenizer = pickle.loads(reader.read()) # Loading tokenizer
        SEQUENCE_SIZE  = 100
        
        # tokenizer
        sequences = tokenizer.texts_to_sequences(df['sequence'])
        x         = pad_sequences(sequences, maxlen=SEQUENCE_SIZE, padding='post', truncating='post', value=0)

        # predict
        predict = model.predict(x) 
        result  = np.argmax(predict,axis=1).tolist()
        return result

    def spgn(df): # gram negative signal peptide function model
        
        model = keras.models.load_model(
            "exe/model/spgn.keras", custom_objects={"TransformerBlock": TransformerBlock}
            ) # Loading model 

        with open('exe/tokenizer/spgn.pkl', 'rb') as reader:
            tokenizer = pickle.loads(reader.read()) # Loading tokenizer
            
        SEQUENCE_SIZE  = 100    

        # tokenizer
        sequences = tokenizer.texts_to_sequences(df['sequence'])
        x         = pad_sequences(sequences, maxlen=SEQUENCE_SIZE, padding='post', truncating='post', value=0)

        # predict
        predict = model.predict(x) 
        predict = predict >= 0.5
        return predict.tolist()


    def spp(df): # gram positive signal pepitede function model
        
        model = keras.models.load_model(
            "exe/model/spp.keras", custom_objects={"TransformerBlock": TransformerBlock}
            ) # Loading model 

        with open('exe/tokenizer/spp.pkl', 'rb') as reader:
            tokenizer = pickle.loads(reader.read()) # Loading tokenizer
        
        SEQUENCE_SIZE  = 100

        # tokenizer
        sequences = tokenizer.texts_to_sequences(df['sequence'])
        x         = pad_sequences(sequences, maxlen=SEQUENCE_SIZE, padding='post', truncating='post', value=0)

        # predict
        predict = model.predict(x) 
        predict = predict >= 0.5
        return predict.tolist()

    def bbn(df): # gram negative beta barrel function model
        
        model = keras.models.load_model(
            "exe/model/bbn.keras", custom_objects={"TransformerBlock": TransformerBlock}
            ) # Loading model 

        with open('exe/tokenizer/bbn.pkl', 'rb') as reader:
            tokenizer = pickle.loads(reader.read()) # Loading tokenizer

        SEQUENCE_SIZE  = 1000

        # tokenizer
        sequences = tokenizer.texts_to_sequences(df['sequence'])
        x         = pad_sequences(sequences, maxlen=SEQUENCE_SIZE, padding='post', truncating='post', value=0)

        # predict
        predict = model.predict(x) 
        predict = predict >= 0.5
        return predict.tolist()

    def ahn(df): # gram negative alpha helix function model
        
        model = keras.models.load_model(
            "exe/model/ahn.h5", custom_objects={"TransformerBlock": TransformerBlock}, compile=False
            ) # Loading model
        model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

        with open('exe/tokenizer/ahn.pkl', 'rb') as reader:
            tokenizer = pickle.loads(reader.read()) # Loading tokenizer

        SEQUENCE_SIZE  = 3000

        # tokenizer
        sequences = tokenizer.texts_to_sequences(df['sequence'])
        x         = pad_sequences(sequences, maxlen=SEQUENCE_SIZE, padding='post', truncating='post', value=0)

        # predict
        predict = model.predict(x) 
        predict = predict >= 0.5
        return predict.tolist()
