from trackc.pl.bed import bed_track
from trackc.pl.bigwig import bw_compartment, bw_track
from trackc.pl.gene import gene_track
from trackc.pl.links import links_track
from trackc.pl.mapc import getData2Map, mapC
from trackc.pl.mapc_markline import mapc_markline
from trackc.pl.scale import multi_scale_track, scale_track
from trackc.pl.vhighlight import vhighlight
from trackc.pl.Virtual4C import virtual4C
from trackc.pl.zoomin import zoomin
from trackc.pl.bedgraph_matrix import bdgmat_track

__doc__ = """\
Plotting API
============

.. currentmodule:: trackc

.. note::
   See the :ref:`API` section for all important plotting configurations.

.. _pl-generic: 

Generic
-------

.. autosummary::
   :toctree: .
   pl.getData2Map
   pl.mapC
   pl.gene_track
   pl.mapc_markline
   pl.virtual4C
   pl.links_track
   pl.bw_track
   pl.bw_compartment
   pl.bed_track
   pl.multi_scale_track
   pl.vhighlight
   pl.bdgmat_track
"""
