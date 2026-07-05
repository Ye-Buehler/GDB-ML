===============
Reproducibility
===============

This document describes how the files and scripts in this repository map to the
workflow described in the manuscript. All paths below are relative to the root
of this repository.


Software Environments
=====================

Transformer
--------------------

The transformer workflow uses OpenNMT-py. Install OpenNMT-py from pip:

.. code-block:: bash

    pip install OpenNMT-py

or from source:

.. code-block:: bash

    git clone https://github.com/OpenNMT/OpenNMT-py.git
    cd OpenNMT-py
    pip install -e .

OpenNMT-py requires Python >= 3.8 and PyTorch >= 2.0, < 2.2. If a
``MemoryError`` occurs during installation, retry with ``--no-cache-dir``.
Optional OpenNMT-py features can be installed with:

.. code-block:: bash

    pip install -r requirements.opt.txt


Generative models
-----------------

The Conda environment for the generative model workflow is provided at:

.. code-block:: text

    generative_models/environment-py39.yaml

Create and activate it with:

.. code-block:: bash

    conda env create -f generative_models/environment-py39.yaml
    conda activate reinvent-gdb13-py39

PySpark requires Java. Please install a JDK, for example JDK 11 or 17, and make
sure ``JAVA_HOME`` is set before running scripts that use PySpark.



Repository File Map
===================

Transformer train and validation files
--------------------------------------

The manuscript uses compact names such as ``src_train.txt`` and
``tgt_train.txt``. In this repository, the corresponding files are split into
parts to keep individual files manageable.

.. list-table::
   :header-rows: 1

   * - Manuscript name
     - Repository files
     - Meaning
   * - ``src_train.txt``
     - ``transformer/gdb20_data/shuffled_train_keys_part_*_canonical_concatenated_tokenized.txt``
     - Tokenized source graph strings for training.
   * - ``tgt_train.txt``
     - ``transformer/gdb20_data/shuffled_train_values_part_*_canonical_concatenated_tokenized.txt``
     - Tokenized target molecule strings for training.
   * - ``src_val.txt``
     - ``transformer/gdb20_data/shuffled_val_keys_part_*_canonical_concatenated_tokenized.txt``
     - Tokenized source graph strings for validation.
   * - ``tgt_val.txt``
     - ``transformer/gdb20_data/shuffled_val_values_part_*_canonical_concatenated_tokenized.txt``
     - Tokenized target molecule strings for validation.

If manuscript-style filenames are desired, concatenate the split files in sorted
part order:

.. code-block:: bash

    cat transformer/gdb20_data/shuffled_train_keys_part_*_canonical_concatenated_tokenized.txt \
        > transformer/gdb20_data/src_train.txt
    cat transformer/gdb20_data/shuffled_train_values_part_*_canonical_concatenated_tokenized.txt \
        > transformer/gdb20_data/tgt_train.txt
    cat transformer/gdb20_data/shuffled_val_keys_part_*_canonical_concatenated_tokenized.txt \
        > transformer/gdb20_data/src_val.txt
    cat transformer/gdb20_data/shuffled_val_values_part_*_canonical_concatenated_tokenized.txt \
        > transformer/gdb20_data/tgt_val.txt


Generative model train and validation files
-------------------------------------------

The generative model input files are located at:

.. code-block:: text

    generative_models/gdb20_data/1M_node1-17_train.txt
    generative_models/gdb20_data/1M_node1-17_validation.txt
    generative_models/gdb20_data/1M_node18_train.txt
    generative_models/gdb20_data/1M_node18_validation.txt
    generative_models/gdb20_data/1M_node19_train.txt
    generative_models/gdb20_data/1M_node19_validation.txt
    generative_models/gdb20_data/1M_node20_train_1.txt
    generative_models/gdb20_data/1M_node20_train_2.txt
    generative_models/gdb20_data/1M_node20_validation_1.txt
    generative_models/gdb20_data/1M_node20_validation_2.txt


Graph Selection Before Transformer Training
===========================================

The graph-selection utilities are implemented in:

.. code-block:: text

    src/gdb_ml/graph_mapping.py

The main class is ``GraphMapping``. The relevant default values are:

.. code-block:: python

    MOLS_PER_GRAPH = 10
    TOTAL_DATAPOINTS = 2000000

The function ``GraphMapping.datapoints_split`` expects a graph summary table as
input. The table must contain at least these columns:

.. code-block:: text

    Key
    Number of Values

The expected meaning is:

* ``Key``: the graph string used as the transformer source/input.
* ``Number of Values``: the number of molecules associated with that graph.

The selection procedure is:

1. Assign ``MOLS_PER_GRAPH`` as the requested number of datapoints per graph.
2. Keep the first ``TOTAL_DATAPOINTS / MOLS_PER_GRAPH`` graph rows from the
   input table. With the defaults above, this corresponds to 200,000 graph rows.
3. Assign rows to train, validation, and test using the repeated pattern
   ``17 train : 2 validation : 1 test``.
4. Sort the train, validation, and test subsets by ``Number of Values`` in
   descending order.
5. For each selected graph, keep up to ``MOLS_PER_GRAPH`` molecules. If a graph
   has fewer than ``MOLS_PER_GRAPH`` associated molecules, all available
   molecules are retained.

The function ``GraphMapping.check_mols_from_graph`` then reads JSON dictionaries
mapping graph keys to molecule lists, merges entries for selected graph keys,
and truncates each selected graph to at most ``MOLS_PER_GRAPH`` molecules.



Transformer Preprocessing
=========================

After preparing the source and target files, run OpenNMT-py preprocessing. The
README gives the full command; the core inputs are:

.. code-block:: bash

    python transformer/preprocess.py \
        -train_src transformer/gdb20_data/src_train.txt \
        -train_tgt transformer/gdb20_data/tgt_train.txt \
        -valid_src transformer/gdb20_data/src_val.txt \
        -valid_tgt transformer/gdb20_data/tgt_val.txt \
        -save_data transformer/data/voc_exp36/Preprocessed \
        -src_seq_length 3000 \
        -tgt_seq_length 3000 \
        -src_vocab_size 3000 \
        -tgt_vocab_size 3000 \
        -share_vocab \
        -lower

Option meanings:

.. list-table::
   :header-rows: 1

   * - Option
     - Purpose
   * - ``-train_src``
     - Tokenized source graph strings used as model inputs during training.
   * - ``-train_tgt``
     - Tokenized target molecule strings used as model outputs during training.
   * - ``-valid_src``
     - Tokenized source graph strings used for validation.
   * - ``-valid_tgt``
     - Tokenized target molecule strings used for validation.
   * - ``-save_data``
     - Prefix/path for the OpenNMT-py preprocessed dataset.
   * - ``-src_seq_length``
     - Maximum source sequence length retained during preprocessing.
   * - ``-tgt_seq_length``
     - Maximum target sequence length retained during preprocessing.
   * - ``-src_vocab_size``
     - Maximum source vocabulary size.
   * - ``-tgt_vocab_size``
     - Maximum target vocabulary size.
   * - ``-share_vocab``
     - Use a shared source and target vocabulary.
   * - ``-lower``
     - Lowercase tokens during preprocessing.


Transformer Training
====================

The README gives the full training command. The most important options are:

.. list-table::
   :header-rows: 1

   * - Option
     - Purpose
   * - ``-data``
     - Prefix/path produced by the preprocessing step.
   * - ``-save_model``
     - Output prefix for saved transformer checkpoints.
   * - ``-seed``
     - Random seed used for reproducibility.
   * - ``-train_steps``
     - Number of optimization steps.
   * - ``-save_checkpoint_steps``
     - Save a checkpoint every N steps.
   * - ``-keep_checkpoint``
     - Number of checkpoints to keep.
   * - ``-batch_size``
     - Batch size. In this workflow it is interpreted with ``-batch_type tokens``.
   * - ``-batch_type tokens``
     - Batch examples by token count rather than by sentence count.
   * - ``-accum_count``
     - Number of gradient accumulation steps.
   * - ``-optim adam``
     - Use the Adam optimizer.
   * - ``-learning_rate``
     - Initial learning-rate scale used with the Noam schedule.
   * - ``-decay_method noam``
     - Use the Noam learning-rate decay schedule.
   * - ``-warmup_steps``
     - Number of warmup steps for the Noam schedule.
   * - ``-layers``
     - Number of encoder and decoder layers.
   * - ``-rnn_size``
     - Hidden size used by the model implementation.
   * - ``-word_vec_size``
     - Token embedding dimension.
   * - ``-encoder_type transformer``
     - Use a transformer encoder.
   * - ``-decoder_type transformer``
     - Use a transformer decoder.
   * - ``-heads``
     - Number of transformer attention heads.
   * - ``-transformer_ff``
     - Feed-forward layer size in transformer blocks.
   * - ``-dropout``
     - Dropout probability.
   * - ``-valid_steps``
     - Run validation every N training steps.
   * - ``-early_stopping``
     - Stop training after this many validations without improvement.
   * - ``-gpu_ranks``
     - GPU rank used for training.
   * - ``-tensorboard``
     - Enable TensorBoard logging.
   * - ``-tensorboard_log_dir``
     - Directory where TensorBoard logs are written.


Transformer Generation
======================

Generation uses ``transformer/translate.py``:

.. code-block:: bash

    python transformer/translate.py \
        -model "$MODEL_PATH" \
        -src "$SRC_FILE" \
        -output "$OUTPUT_FILE" \
        -batch_size 64 \
        -replace_unk \
        -max_length 1000 \
        -log_probs \
        -beam_size 300 \
        -n_best 300

Option meanings:

.. list-table::
   :header-rows: 1

   * - Option
     - Purpose
   * - ``-model``
     - Trained transformer checkpoint.
   * - ``-src``
     - Source graph file used for generation.
   * - ``-output``
     - File where generated target strings are written.
   * - ``-batch_size``
     - Number of source examples processed per batch.
   * - ``-replace_unk``
     - Replace unknown tokens where possible.
   * - ``-max_length``
     - Maximum generated sequence length.
   * - ``-log_probs``
     - Write log probabilities for generated sequences.
   * - ``-beam_size``
     - Number of hypotheses maintained during beam search.
   * - ``-n_best``
     - Number of generated hypotheses retained per source example.


Generative Model Workflow
=========================

The generative model scripts are located in ``generative_models/``:

.. code-block:: text

    create_randomized_smiles.py
    create_model.py
    train_model.py
    sample_from_model.py
    calculate_nlls.py

Example workflow:

.. code-block:: bash

    cd generative_models
    mkdir -p node18_randomized/models

    ./create_randomized_smiles.py \
        -i gdb20_data/1M_node18_train.txt \
        -o node18_randomized/training \
        -n 100

    ./create_randomized_smiles.py \
        -i gdb20_data/1M_node18_validation.txt \
        -o node18_randomized/validation \
        -n 100

    ./create_model.py \
        -i node18_randomized/training/001.smi \
        -o node18_randomized/models/model.empty

    ./train_model.py \
        -i node18_randomized/models/model.empty \
        -o node18_randomized/models/model.trained \
        -s node18_randomized/training \
        -e 100 \
        --lrm ada \
        --csl node18_randomized/tensorboard \
        --csv node18_randomized/validation \
        --csn 75000

    ./sample_from_model.py \
        -m node18_randomized/models/model.trained.100 \
        -n 1000000 \
        --with-nll \
        -o output.txt

Option meanings for ``create_randomized_smiles.py``:

.. list-table::
   :header-rows: 1

   * - Option
     - Purpose
   * - ``-i``, ``--input-smi-path``
     - Input SMILES file to randomize.
   * - ``-o``, ``--output-smi-folder-path``
     - Output folder where randomized SMILES files are written.
   * - ``-n``, ``--num-files``
     - Number of randomized SMILES files to create. Files are numbered from
       ``000.smi``.
   * - ``-r``, ``--random-type``
     - Randomization mode. Supported values are ``restricted`` and
       ``unrestricted``.
   * - ``-s``, ``--smiles-type``
     - Input/output SMILES representation type, for example ``smiles`` or a
       supported DeepSMILES mode.
   * - ``-p``, ``--num-partitions``
     - Number of Spark partitions used when processing the input file.

Option meanings for ``create_model.py``:

.. list-table::
   :header-rows: 1

   * - Option
     - Purpose
   * - ``-i``, ``--input-smiles-path``
     - SMILES file used to build the model vocabulary.
   * - ``-o``, ``--output-model-path``
     - Path/prefix where the empty initialized model is saved.
   * - ``-l``, ``--num-layers``
     - Number of recurrent neural-network layers.
   * - ``-s``, ``--layer-size``
     - Hidden size of each recurrent layer.
   * - ``-e``, ``--embedding-layer-size``
     - Embedding-layer dimension.
   * - ``-d``, ``--dropout``
     - Dropout applied between GRU layers.
   * - ``--ln``, ``--layer-normalization``
     - Enable layer normalization on GRU output.
   * - ``--max-sequence-length``
     - Maximum generated/encoded sequence length.

Option meanings for ``train_model.py``:

.. list-table::
   :header-rows: 1

   * - Option
     - Purpose
   * - ``-i``, ``--input-model-path``
     - Input model file, usually the empty model created by
       ``create_model.py``.
   * - ``-o``, ``--output-model-prefix-path``
     - Output model prefix. The epoch number is appended when checkpoints are
       saved.
   * - ``-s``, ``--training-set-path``
     - Training SMILES file or directory containing multiple ``.smi`` files.
   * - ``-e``, ``--epochs``
     - Number of training epochs.
   * - ``-b``, ``--batch-size``
     - Number of molecules processed per batch.
   * - ``--sen``, ``--save-every-n-epochs``
     - Save the model every N epochs.
   * - ``--clip-gradients``
     - Clip gradients to the given norm.
   * - ``--lrm``, ``--learning-rate-mode``
     - Learning-rate schedule mode. Supported values are ``exp`` and ``ada``.
   * - ``--lrs``, ``--learning-rate-start``
     - Starting learning rate.
   * - ``--lrmin``, ``--learning-rate-min``
     - Minimum learning rate; training stops when this value is reached.
   * - ``--lrg``, ``--learning-rate-gamma``
     - Multiplicative factor used when lowering the learning rate.
   * - ``--lrt``, ``--learning-rate-step``
     - Number of epochs between learning-rate changes in exponential mode.
   * - ``--lrth``, ``--learning-rate-threshold``
     - Threshold used to lower the learning rate in adaptive mode.
   * - ``--lras``, ``--learning-rate-average-steps``
     - Number of previous metric values used for adaptive learning-rate
       averaging.
   * - ``--lrp``, ``--learning-rate-patience``
     - Number of unimproved steps before lowering the learning rate in adaptive
       mode.
   * - ``--csf``, ``--collect-stats-frequency``
     - Collect validation/statistics every N epochs.
   * - ``--csl``, ``--collect-stats-log-path``
     - TensorBoard/statistics output directory.
   * - ``--csv``, ``--collect-stats-validation-set-path``
     - Validation SMILES file or directory used when collecting statistics.
   * - ``--csn``, ``--collect-stats-sample-size``
     - Number of SMILES sampled from the model when collecting statistics.
   * - ``--csw``, ``--collect-stats-with-weights``
     - Store model weight matrices when collecting statistics.
   * - ``--csst``, ``--collect-stats-smiles-type``
     - SMILES representation used for statistics, for example ``smiles`` or a
       supported DeepSMILES mode.

Option meanings for ``sample_from_model.py``:

.. list-table::
   :header-rows: 1

   * - Option
     - Purpose
   * - ``-m``, ``--model-path``
     - Trained model checkpoint used for sampling.
   * - ``-o``, ``--output-smiles-path``
     - Output file. If omitted, samples are written to standard output.
   * - ``-n``, ``--num``
     - Number of SMILES strings to sample.
   * - ``--with-nll``
     - Write the negative log likelihood in a second column after each SMILES.
   * - ``-b``, ``--batch-size``
     - Sampling batch size.
   * - ``--use-gzip``
     - Compress the output file with gzip.


Reproducibility Boundary
========================

This repository contains the source code, scripts, released tokenized
transformer files, generative-model input files, and trained model artifacts
needed to rerun the published workflows from the released intermediate data.
Regenerating every intermediate object before graph selection additionally
requires the graph summary table and graph-to-molecule dictionaries described
above. These inputs should be made available with the repository or with the
linked Zenodo data release for full regeneration from the earliest graph
selection stage.
