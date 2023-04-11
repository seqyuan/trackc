.. image:: https://zenodo.org/badge/572302042.svg
   :target: https://zenodo.org/badge/latestdoi/572302042

trackc

install
=======

.. code-block:: bash

    pip install trackc

core usage
===========

.. code-block:: python

	import trackc

	kegg_metabolism = sigc.metabolism_sigs(resources='KEGG')
	# other custom signature gene sets also support

	df = pd.read_table("cells_X_genes.mat", header=0, index_col=0)

	sig_mtx = sigc.sigc_score(df, kegg_metabolism, method="AUCell")
	# sig_mtx: cells_X_signatures

signature example
==================

============== ============ =======
name           description  member
============== ============ =======
Glycolysis     00010        HK3
Glycolysis     00010        HK1
Glycolysis2    describ2     geneX 
============== ============ =======

Documentation
==================

Extensive documentation and tutorials are available at https://trackc.readthedocs.io

