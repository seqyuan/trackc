.. _tutorials:

**********
Tutorials
**********

In this section we present several tutorials on using sigc.

Basic sigc analysis
======================

**Running sigc step by step**

***core code*** 

import sigc

kegg_metabolism = sigc.metabolism_sigs(resources='KEGG')

df = pd.read_table("cells_X_genes.mat", header=0, index_col=0)

sig_mtx = sigc.sigc_score(df, kegg_metabolism, method="AUCell")


.. toctree::
    :hidden:
    :maxdepth: 2
    
    sigc_main.ipynb
    seurat2scanpy.ipynb
