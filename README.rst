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
