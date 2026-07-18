======
Sampling a Trillion-Sized GDB-20 Database of Drug-Like Molecules by Generative Artificial Intelligence
======


Thank you for your interest in this repository, which complements the publication 
"`Sampling a Trillion-Sized GDB-20 Database of Drug-Like Molecules by Generative Artificial Intelligence <https://chemrxiv.org/doi/full/10.26434/chemrxiv.15000288/v1>`_".

For a step-by-step map from the manuscript workflow to the repository files,
software environments, graph-selection procedure, and executable scripts, see
`REPRODUCIBILITY.rst <REPRODUCIBILITY.rst>`_.

.. image:: docs/GA.jpg
   :alt: GA
   :align: center
   :width: 400px

GDB-20s Database:
========================================================================================

**Zenodo Part 1:**

.. code-block:: bash

   # See https://zenodo.org/records/17368725.

**Zenodo Part 2:**

.. code-block:: bash

    # See https://zenodo.org/records/17375415.

Main Codes and Files for GDB-20 Machine-Learning-Based Project:
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

General Usage:
========================================================================================

Transformer Examples:
-----------------------

**(1) Create a Conda environment with the provided `.yaml` file and activate it:**

.. code-block:: bash

    conda env create -f environment-gdb20.yaml
    conda activate gdb20

**(2) Follow the pipeline and tokenize the SMILES:**

.. code-block:: bash
   
   # See transformer/pipeline.ipynb
   
   
**(3) OpenNMT-py Installation:**

.. code-block:: bash

    conda env create -f transformer/environment-opennmt.yaml
    conda activate opennmt
    git clone https://github.com/reymond-group/OpenNMT-py.git
    cd OpenNMT-py
    git checkout Enzymatic_Transformer
    pip install -e .
    cd ..

**(4) Preprocess the data:**

.. code-block:: bash

    # Concatenate split files into manuscript-style filenames
    cat transformer/gdb20_data/shuffled_train_keys_part_*_canonical_concatenated_tokenized.txt \
        > transformer/gdb20_data/src_train.txt
    cat transformer/gdb20_data/shuffled_train_values_part_*_canonical_concatenated_tokenized.txt \
        > transformer/gdb20_data/tgt_train.txt
    cat transformer/gdb20_data/shuffled_val_keys_part_*_canonical_concatenated_tokenized.txt \
        > transformer/gdb20_data/src_val.txt
    cat transformer/gdb20_data/shuffled_val_values_part_*_canonical_concatenated_tokenized.txt \
        > transformer/gdb20_data/tgt_val.txt

    cat \
    "transformer/gdb20_model/test36_model_step_55000_part_aa.pt" \
    "transformer/gdb20_model/test36_model_step_55000_part_ab.pt" \
    "transformer/gdb20_model/test36_model_step_55000_part_ac.pt" \
    > "transformer/gdb20_model/test36_model_step_55000.full.pt"

    # Define variables
    data_dir="gdb20_data"
    dataset="test36"
    experiment="exp36"

    batchsize=6144
    dropout=0.1
    rnnsize=384
    wordvecsize=384
    learnrate=2
    layers=4
    heads=8

    MODEL_PATH="gdb20_model/test36_model_step_55000.full.pt"
    SRC_FILE="gdb20_data/src_val.txt"
    OUTPUT_FILE="experiments/test36_model_step_55000_predictions.txt"

    # Run the remaining transformer commands from this directory.
    cd transformer

    mkdir -p data/voc_${experiment}

    # Run preprocessing
    python preprocess.py \
        -train_src "${data_dir}/src_train.txt" \
        -train_tgt "${data_dir}/tgt_train.txt" \
        -valid_src "${data_dir}/src_val.txt" \
        -valid_tgt "${data_dir}/tgt_val.txt" \
        -save_data data/voc_${experiment}/Preprocessed \
        -src_seq_length 3000 -tgt_seq_length 3000 \
        -src_vocab_size 3000 -tgt_vocab_size 3000 \
        -share_vocab -lower


**(5) Train the transformer model:**

.. code-block:: bash

    # Remove the line "-gpu_ranks 0 \"" when training on CPU
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


**(6) Molecular generation:**

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


Generative Models Examples:
-----------------------------

System Dependency:
~~~~~~~~~~~~~~~~~~

PySpark requires Java. Please install a JDK, e.g. JDK 11 or 17, and make sure
``JAVA_HOME`` is set before running scripts that use PySpark.



**(1) Activate the Conda environment:**

.. code-block:: bash

    conda activate gdb20

    # While gdb20 is active, install Java:
    conda install -c conda-forge openjdk=17

**(2) Create a working directory:**

.. code-block:: bash

    # Run the remaining RNN commands from this directory.
    cd generative_models
    
    mkdir -p node18_randomized/models

**(3) Create random SMILES:**

.. code-block:: bash

    ./create_randomized_smiles.py -i gdb20_data/1M_node18_train.txt -o node18_randomized/training -n 100
    ./create_randomized_smiles.py -i gdb20_data/1M_node18_validation.txt -o node18_randomized/validation -n 100

**(4) Create a blank model file:**

.. code-block:: bash

    ./create_model.py -i node18_randomized/training/001.smi -o node18_randomized/models/model.empty

Note:
~~~~~~~~~~~~~~~~~~

Steps 5 and 6 require an NVIDIA GPU with CUDA. If using a GPU, replace
``pytorch=2.7.0`` with ``pytorch-gpu=2.7.0`` in ``environment-gdb20.yaml`` before creating the environment.

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

    # To use the bundled pretrained model instead, replace the -m path below with:
    # gdb20_models/model.trained.node18
    ./sample_from_model.py \
        -m node18_randomized/models/model.trained.100 \
        -n 1000000 \
        --with-nll \
        -o output.txt



Original OpenNMT-py:
--------

* If you reuse this code please also cite the underlying code frameworks: "`OpenNMT technical report <https://www.aclweb.org/anthology/P17-4012>`_" and "`Enzymatic_Transformer <https://github.com/reymond-group/OpenNMT-py>`_".

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
