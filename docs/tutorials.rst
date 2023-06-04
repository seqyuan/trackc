##########
Tutorials
##########

GridSpec
========

.. code-block:: python

    import trackc as tc
    fig, axs = tc.make_spec(width=6, height=4, height_ratios=[1,1.5,2], hspace=0.3)

.. image-sg:: /images/tutorials/grid_spec_01.png
   :alt: Neighborhood enrichment
   :srcset: /images/tutorials/grid_spec_01.png
   :class: sphx-glr-single-img


Analysis of Hi-C using :mod:`trackc`
-------------------------------------------------
This section contains tutorials showcasing core trackc functionalities by applying them
to a diverse set of different Hi-C datasets.

input file formats
==================

.. toctree::
   :maxdepth: 2
   :caption: File formats:

   tutorials/fileformats/bed12
   tutorials/fileformats/mcool

Available track types
======================
This section contains various quick tutorials showcasing omics data visualization with :mod:`trackc`.

.. toctree::
   :maxdepth: 2
   :caption: Available track types:

   tutorials/scalebar
   tutorials/contactmap
   tutorials/Virtual4C
   tutorials/gene
   tutorials/bigwig
   tutorials/bed
   tutorials/links
   tutorials/zoomin
   