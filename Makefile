setup:
	conda env create --file environment.yml || conda env update --file environment.yml

training_sp_gram_positive:
	python exe/model_training_signal_peptide_gram_positive.py

training_sp_gram_negative:
	python exe/model_training_signal_peptide_gram_negative.py