.. _Installation:

=====================
Installation
=====================

Clone repository
=================
Clone the github repository

.. code-block:: bash

  git clone https://github.com/Ahituv-lab/Asymmetron


Set up the conda environment
============================

Create a conda virtual enviroment with the necesary dependencies

.. code-block:: bash

  conda config --add channels defaults

  conda config --add channels bioconda

  conda config --add channels conda-forge

  conda create --name asymmetron pybedtools python=3.5 seaborn numpy


Initiate Asymmetron
====================
Activate virtual enviroment

.. code-block:: bash

  conda activate asymmetron
