from trackc.plotting.scale import multi_scale_track, scale_track
from trackc.plotting.bed import bed_track
from trackc.plotting.bigwig import bw_track, bw_compartment
from trackc.plotting.links import links_track
from trackc.plotting.Virtual4C import virtual4C
from trackc.plotting.mapc_markline import mapc_markline
from trackc.plotting.gene import gene_track
from trackc.plotting.mapc import (
    mapC,
    getData2Map
)
from trackc.plotting.zoomin import zoomin

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
"""
