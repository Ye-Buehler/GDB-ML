======
GDB-20 XXX
======


.. image:: https://img.shields.io/pypi/v/gdb_ml.svg
        :target: https://pypi.python.org/pypi/gdb_ml

.. image:: https://readthedocs.org/projects/gdb-ml/badge/?version=latest
        :target: https://gdb-ml.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

Thank you for your interest in this repository, which complements the publication 
"`GDB-20: XXX <https://XXX>`_".

.. image:: https://github.com/Ye-Buehler/XXX.jpg
   :alt: GA
   :align: center
   :width: 400px


Main codes and files for GDB-20 machine-learning-based project:
========================================================================================

.. code-block:: text

    GDB-ML/
    ├── src/
    |   └── gdb_ml/
    |      ├── chem_utils.py
    |      ├── data_processor.py
    |      ├── graph_mapping.py
    |      └── properties_calculator.py
    ├── transformer/
    |      ├── pipeline.ipynb
    |      ├── preprocess.py
    |      ├── train.py
    |      ├── translate.py
    |      ├── gdb20_data
    |      └── gdb20_model
    └── generative_models/
           ├── create_randomized_smiles.py
           ├── create_model.py
           ├── train_model.py
           ├── sample_from_model.py
           ├── calculate_nlls.py
           ├── gdb20_data
           └── gdb20_models

General usage:
========================================================================================

Transformer Examples
-----------------------

**(1) Follow the pipeline and tokenize the SMILES:**

.. code-block:: bash

   # See pipeline.ipynb

**(2) Preprocess the data:**

.. code-block:: bash

    # Activate environment
    conda activate opennmt

    # Define variables
    dataset="test36"
    experiment="exp36"

    batchsize=6144
    dropout=0.1
    rnnsize=384
    wordvecsize=384
    learnrate=2
    layers=4
    heads=8

    mkdir -p data/voc_${experiment}

    # Run preprocessing
    python preprocess.py \
        -train_src data/${dataset}/shuffled_train_keys_canonical_concatenated_tokenized.txt \
        -train_tgt data/${dataset}/shuffled_train_values_canonical_concatenated_tokenized.txt \
        -valid_src data/${dataset}/shuffled_val_keys_canonical_concatenated_tokenized.txt \
        -valid_tgt data/${dataset}/shuffled_val_values_canonical_concatenated_tokenized.txt \
        -save_data data/voc_${experiment}/Preprocessed \
        -src_seq_length 3000 -tgt_seq_length 3000 \
        -src_vocab_size 3000 -tgt_vocab_size 3000 \
        -share_vocab -lower


**(3) Train the Transformer model:**

.. code-block:: bash

    python train.py \
        -data data/voc_${experiment}/Preprocessed \
        -save_model experiments/checkpoints/${experiment}/${dataset}_model \
        -seed 42 \
        -save_checkpoint_steps 500 \
        -keep_checkpoint 50 \
        -train_steps 500000 \
        -param_init 0 \
        -param_init_glorot \
        -max_generator_batches 32 \
        -batch_size ${batchsize} \
        -batch_type tokens \
        -normalization tokens \
        -max_grad_norm 0 \
        -accum_count 4 \
        -optim adam \
        -adam_beta1 0.9 \
        -adam_beta2 0.998 \
        -decay_method noam \
        -warmup_steps 8000 \
        -learning_rate ${learnrate} \
        -label_smoothing 0.0 \
        -layers 4 \
        -rnn_size ${rnnsize} \
        -word_vec_size ${wordvecsize} \
        -encoder_type transformer \
        -decoder_type transformer \
        -dropout ${dropout} \
        -position_encoding \
        -global_attention general \
        -global_attention_function softmax \
        -self_attn_type scaled-dot \
        -heads 8 \
        -transformer_ff 2048 \
        -valid_steps 500 \
        -valid_batch_size 4 \
        -report_every 500 \
        -log_file data/Training_LOG_${experiment}.txt \
        -early_stopping 10 \
        -early_stopping_criteria accuracy \
        -world_size 1 \
        -gpu_ranks 0 \
        -tensorboard \
        -tensorboard_log_dir experiments/Tensorboard/${experiment}/


**(4) Molecular generation:**

.. code-block:: bash

    python translate.py \
        -model "$MODEL_PATH" \
        -src "$SRC_FILE" \
        -output "$OUTPUT_FILE" \
        -batch_size 64 \
        -replace_unk \
        -max_length 1000 \
        -log_probs \
        -beam_size 300 \
        -n_best 300


Generative Models Examples
-----------------------------

**(1) Create a Conda environment with the provided `.yaml` file and activate it:**

.. code-block:: bash

    conda env create -f environment-py39.yaml
    conda activate reinvent-gdb13-py39

**(2) Create a working directory:**

.. code-block:: bash

    mkdir -p node18_randomized/models

**(3) Create random SMILES:**

.. code-block:: bash

    ./create_randomized_smiles.py -i training_sets/1M_node18_train.txt -o node18_randomized/training -n 100
    ./create_randomized_smiles.py -i training_sets/1M_node18_validation.txt -o node18_randomized/validation -n 100

**(4) Create a blank model file:**

.. code-block:: bash

    ./create_model.py -i node18_randomized/training/001.smi -o node18_randomized/models/model.empty

**(5) Train the generative model with specified parameters:**

.. code-block:: bash

    ./train_model.py \
        -i node18_randomized/models/model.empty \
        -o node18_randomized/models/model.trained \
        -s node18_randomized/training \
        -e 100 --lrm ada \
        --csl node18_randomized/tensorboard \
        --csv node18_randomized/validation \
        --csn 75000

**(6) Sample an already trained model for a given number of SMILES (also retrieves log-likelihoods):**

.. code-block:: bash

    ./sample_from_model.py \
        -m node18_randomized/models/model.trained.100 \
        -n 1000000 \
        --with-nll \
        -o output.txt



Original OpenNMT-py:
--------

* If you reuse this code please also cite the underlying code framework: "`OpenNMT technical report <https://www.aclweb.org/anthology/P17-4012>`_" and "`Enzymatic_Transformer <https://github.com/reymond-group/OpenNMT-py>`_".

Original Reinvent-Randomized:
--------

* If you reuse this code please also cite the underlying code framework: "`reinvent-randomized <https://github.com/undeadpixel/reinvent-randomized>`_".

License
--------

* Free software: MIT license


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
