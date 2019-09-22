.. |date| date::

*******************************
PPTC TMB Plots
*******************************

:authors: Jo Lynne Rokita
:contact: Jo Lynne Rokita (rokita@email.chop.edu)
:organization: CHOP
:status: In-process
:date: |date|

.. meta::
   :keywords: pdx, mouse, WES, TMB, 2019
   :description: code to create TMB plots

Introduction
============

Here, we provide scripts to enable reproducible generation of Manuscript Figures 2B and S3B: TMB plots for diagnosis/relapse samples across the cohort, within neuroblastoma, BCP-ALL, T-ALL, and osteosarcoma, and between paired diagnosis/relapse samples. This repo contains code for:

1. Creating TMB plot (Figure S3B) and TMB violin plots (Figure 3B)
2. Creating summary statistic tables included in Table S4: Ages and TMB by histology

Details
=======

- TMB.R
- demog-color-function.R                  
- install-packages.R
- theme.R



Software Requirements
=====================

R 3.4.3

Pipeline
========

.. code-block:: bash

         # How to run:
         # Download github repository in your home directory (~/)
         git clone https://github.com/marislab/pptc-pdx-tmb.git
   
         # Run script to create pie chart
         Rscript ~/pptc-pdx-tmb/R/TMB.R 
         
Results
========

Results will appear in the results folder

