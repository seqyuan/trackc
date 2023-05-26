.. image:: https://zenodo.org/badge/572302042.svg
   :target: https://zenodo.org/badge/latestdoi/572302042


=======
trackc
=======

Track view of chromosome conformation and multi-omics data
===========================================================
**trackc** is a python package to flexible visualization of 3D genome and multi-omics data.
It builds on top of `matplotlib`, from which it allow for flexible adjustments.

Trackc's key applications
--------------------------
- Mark the abnormal interaction formed by the structural variation of the genome.
- Show the 3Dgenome interaction and multi-omics data after rearrangement.
- Flexible and convenient layout for multi-track 

Available track types
---------------------
- multi Hi-C heatmap
- Virtual4C
- bigwig signal
- bed position, ChIP-seq or ATAC-seq peaks, TAD domains, SVs
- mark region on heatmap
- zoomin
- gene
- links, interaction loops
- scalebar

Here is an example of cite du_dynamic_2022 article diagram implemented by trackc.


Installation
=======

.. code-block:: bash

    pip install trackc

core usage
===========

.. code-block:: python

	import trackc as tc

Documentation
=============
Our documentation provide the full list of possible track types and gallary guidelines for users.

Extensive documentation and tutorials are available at https://trackc.readthedocs.io


Citation
========
If you use trackc in your analysis, Please cite trackc as follows:

Zan Yuan, Guoliang Li, Yang Chen.
**trackc: a Package for Flexible Visualization of rearrangement 3D Genome and Multi-omics**, Bioinformatics, Volume 48, Issue W1, 1 July 2023, Pages W177–W184, https://doi.org/10.1093/nar/gkaa220
