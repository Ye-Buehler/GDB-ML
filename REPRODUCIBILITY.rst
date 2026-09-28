===============
Reproducibility
===============

This document describes how the files and scripts in this repository map to the
workflow described in the manuscript. All paths below are relative to the root
of this repository. We do not claim that the dependency versions specified in this repository exactly reconstruct the original study environment; they define a precise environment for future reproductions of the released workflow.

End-to-End Execution Sequence
=====================

The repository components are connected in the following order:

1. Create and activate the required Conda environment.

2. Install the bundled OpenNMT-py implementation from ``transformer/onmt/``.

3. Concatenate the split files in ``transformer/gdb20_data/`` to form
   ``src_train.txt``, ``tgt_train.txt``, ``src_val.txt`` and
   ``tgt_val.txt``.

4. Run ``transformer/preprocess.py`` with the settings given above to
   produce the OpenNMT preprocessed training and validation data.

5. Run ``transformer/train.py`` or load the released checkpoint from
   ``transformer/gdb20_model/``, and then use ``transformer/translate.py``
   with the specified beam-search settings to generate molecular SMILES
   from the source graphs.

6. Detokenize the transformer output following ``transformer/pipeline.ipynb``. 

7. Retain RDKit-valid canonical SMILES, restrict the molecules to the target heavy-atom-count range and remove duplicate structures.

8. For the RNN workflow, use the five HAC-stratified training and validation
   datasets supplied in ``generative_models/gdb20_data/``. Randomize their SMILES with
   ``create_randomized_smiles.py``, initialize and train the models with
   ``create_model.py`` and ``train_model.py``, or use the released
   checkpoints in ``generative_models/gdb20_models/``.

9. Generate RNN SMILES with ``sample_from_model.py`` and apply the same validity, canonicalization, heavy-atom-count and within-model deduplication procedures.

10. Combine the unique transformer and RNN outputs and remove structures occurring in both outputs. The resulting union constitutes GDB-20s.

Note that this procedure can generate a new molecular collection using the same methodology but does not exactly regenerate the released 12-billion-molecule GDB-20s collection. The exact released GDB-20s collection is provided separately through the Zenodo records linked in ``README.rst``.


Software Environments
=====================

Transformer pipeline and generative models
------------------------------

``environment-gdb20.yaml`` is the canonical environment for the pipeline
notebook, the ``gdb_ml`` helpers, and the generative-model workflow. It uses
Python 3.10 and pins the principal scientific and machine-learning packages
used for the released workflow.

Create and activate the environment from the repository root:

.. code-block:: bash

    conda env create -f environment-gdb20.yaml
    conda activate gdb20

``create_randomized_smiles.py`` uses PySpark and therefore requires Java.
Install a JDK after activating ``gdb20`` and verify it before running that
script:

.. code-block:: bash

    conda install -c conda-forge openjdk=17
    java -version

The YAML uses Conda-forge builds of PyTorch so that PyTorch and the Conda
scientific stack share the same OpenMP runtime. On Linux, Conda can select a
CUDA build when a compatible NVIDIA driver is detected; on macOS and
CPU-only systems it selects a CPU build. Check the result with:

.. code-block:: bash

    python -c "import torch; print(torch.__version__); print(torch.cuda.is_available())"

The released generative training and sampling code contains explicit
``.cuda()`` calls and therefore requires a CUDA-enabled PyTorch build unless
those calls are replaced by device-aware CPU/MPS handling.


Transformer training and generation
-----------

The transformer workflow must be kept separate from the gdb20 environment 
because it requires an older PyTorch stack. 
Create the pinned environment supplied with this repository, 
and install the included OpenNMT-py implementation:

.. code-block:: bash

    conda env create -f transformer/environment-opennmt.yaml
    conda activate opennmt
    git clone https://github.com/reymond-group/GDB-ML.git
    cd GDB-ML
    git checkout v1.0.1
    pip install -e ./transformer
    cd ..

The repository's transformer commands and option names correspond to 
the bundled OpenNMT-py implementation in transformer/onmt 
and may differ from current OpenNMT-py releases. 
On a CPU-only machine, omit -gpu_ranks 0 from the training command in README.rst.

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


Transformer Preprocessing
=========================

After preparing the source and target files, activate the separate OpenNMT
environment and run preprocessing from ``transformer/``. The README gives the
full workflow; the core inputs are:

.. code-block:: bash

    conda activate opennmt
    cd transformer
    mkdir -p data/voc_exp36

    python preprocess.py \
        -train_src gdb20_data/src_train.txt \
        -train_tgt gdb20_data/tgt_train.txt \
        -valid_src gdb20_data/src_val.txt \
        -valid_tgt gdb20_data/tgt_val.txt \
        -save_data data/voc_exp36/Preprocessed \
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

Run ``train.py`` from ``transformer/`` in the ``opennmt`` environment. The
README gives the full training command. Its ``-gpu_ranks 0`` option selects the
first CUDA GPU and must be omitted for CPU training. The most important options
are:

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

The released transformer checkpoint is split into three parts. Reconstruct it
from the repository root before generation:

.. code-block:: bash

    cat \
        transformer/gdb20_model/test36_model_step_55000_part_aa.pt \
        transformer/gdb20_model/test36_model_step_55000_part_ab.pt \
        transformer/gdb20_model/test36_model_step_55000_part_ac.pt \
        > transformer/gdb20_model/test36_model_step_55000.full.pt

Generation then uses ``translate.py`` from ``transformer/``:

.. code-block:: bash

    cd transformer

    MODEL_PATH="gdb20_model/test36_model_step_55000.full.pt"
    SRC_FILE="gdb20_data/src_val.txt"
    OUTPUT_FILE="experiments/test36_model_step_55000_predictions.txt"

    mkdir -p experiments

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

    conda activate gdb20
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
        -i node18_randomized/training/000.smi \
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

``create_randomized_smiles.py`` numbers its outputs starting at ``000.smi``.
With the output prefix and epoch count above, the final newly trained checkpoint
is ``node18_randomized/models/model.trained.100``. To sample the released node-18
checkpoint instead, use:

.. code-block:: bash

    ./sample_from_model.py \
        -m gdb20_models/model.trained.node18 \
        -n 1000000 \
        --with-nll \
        -o output.txt

The training and sampling scripts use explicit CUDA tensor placement. Confirm
that ``torch.cuda.is_available()`` is ``True`` before running them unchanged.
The data-randomization and blank-model creation steps do not require CUDA.

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
     - Dropout applied to embeddings and between LSTM layers.
   * - ``--ln``, ``--layer-normalization``
     - Legacy command-line option. The current ``create_model.py`` does not
       pass this value to the LSTM implementation.
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

This repository supports rerunning the published training and sampling
workflows from the released intermediate files. It includes:

* source code and executable scripts;
* tokenized transformer input files in ``transformer/gdb20_data/``;
* generative-model input files in ``generative_models/gdb20_data/``;
* trained model artifacts in ``transformer/gdb20_model/`` and
  ``generative_models/gdb20_models/``.

The repository does not contain the earliest intermediate files required to reproduce the graph-selection step. The detailed procedure is described in the manuscript; the following section summarizes the relevant workflow:

* Generate planar molecular graphs with up to 20 nodes using GENG, excluding three- and four-membered rings.
* Retain graphs satisfying the reported structural criteria: at most three rings, no node shared by three rings, at most one seven- or eight-membered ring, no larger rings, and at least 40% divalent nodes (`MC1 < 0.6`). These are referred to as GDB-20 graphs.
* From GDB-11, GDB-13, and GDB-17, extract graphs from molecules satisfying the polarity and functional-group criteria.
* Group molecules by graph, rank the graphs by frequency, and retain at most 300 molecules per graph. 
* Split the graph groups—not individual molecules—into 80% training and 20% validation sets. No graph category is shared between the two sets.
* Concatenate the SMILES of pairs of molecules for two thirds of the hydrocarbon training examples. Apply the same concatenation to their corresponding molecular SMILES, so that each input remains aligned with its target output.
* Use the selected GDB-20 graphs, up to 20 nodes, as the unlabeled generation set.
