from . import palettes

from ._utils import matrix
from ._utils import timeseries, timeseries_subplot, timeseries_as_heatmap


__doc__ = """\
Plotting API
============

.. currentmodule:: trackc

.. note::
   See the :ref:`settings` section for all important plotting configurations.

.. _pl-generic: 

Generic
-------

.. autosummary::
   :toctree: .

   pl.mapC
   pl.trackBw
   pl.trackBed
   pl.trackGene
   pl.arc
   pl.insulation
   pl.compartment



Classes
-------

These classes allow fine tuning of visual parameters. 

.. autosummary::
   :toctree: .

    pl.DotPlot
    pl.MatrixPlot
    pl.StackedViolin


Tools
-----

Methods that extract and visualize tool-specific annotation in an
:class:`~anndata.AnnData` object.  For any method in module ``tl``, there is
a method with the same name in ``pl``.

PCA
~~~
.. autosummary::
   :toctree: .

   pl.pca
   pl.pca_loadings
   pl.pca_variance_ratio
   pl.pca_overview

Embeddings
~~~~~~~~~~
.. autosummary::
   :toctree: .

   pl.tsne
   pl.umap
   pl.diffmap
   pl.draw_graph
   pl.spatial
   pl.embedding

Compute densities on embeddings.

.. autosummary::
   :toctree: .

   pl.embedding_density

Branching trajectories and pseudotime, clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Visualize clusters using one of the embedding methods passing ``color='louvain'``.

.. autosummary::
   :toctree: .

   pl.dpt_groups_pseudotime
   pl.dpt_timeseries
   pl.paga
   pl.paga_path
   pl.paga_compare


Simulations
~~~~~~~~~~~
.. autosummary::
   :toctree: .

   pl.sim
"""
