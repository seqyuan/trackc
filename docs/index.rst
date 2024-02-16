trackc - track view of chromosome conformation and multi-omics data
===================================================================

**trackc** is a Python Package to Flexible Visualization of 3D Genome and Multiomics.
It is built based on `matplotlib`, from which it allows for flexible adjustments.

trackc's key applications
--------------------------
- Mark the abnormal interaction formed by the structural variation of the genome.
- Show the three-dimensional genome interaction and multi-omics data after rearrangement.
- Flexible and convenient layout for multi-track 

Getting started with trackc
----------------------------
- Browse :doc:`Available track-types </track_types/index>` and :doc:`examples </cli/index>`.


Usage modes
-----------
trackc provides two usage modes: command line and API.

Command Line:
^^^^^^^^^^^^^
Users can utilize TrackC via the command line interface (CLI). Upon installation, users can directly invoke TrackC commands to access its functionality. The CLI can be configured using a configuration file. In addition to its core features, TrackC integrates common Hi-C data processing commands.

API:
^^^^
Users can interact with TrackC programmatically through its API. This allows for seamless integration of TrackC's functionalities into custom software applications.

To explore the available functionalities and commands, users can refer to the documentation or utilize the trackc --help command.

Below are the currently integrated functional modules:

- trackc cli
- trackc gtf2bed
- trackc tadscore
- trackc fa2GC
- trackc compartment
- trackc insulationscore


:doc:`Hi-C data analysis guides </analysis_guide/index>` Common Hi-C data analysis tutorials provide input files that can be used for plotting with TrackC.


Contributing to trackc
-----------------------
We are happy about any contributions! Welcome to talk at github issue.

.. toctree::
   :caption: General:
   :maxdepth: 1

   install
   track_types/index
   gallery/index
   api/index
   cli/index
   analysis_guide/index
   developer/index

.. _github: https://github.com/seqyuan/trackc


Citation
---------

Please cite trackc as follows:

**trackc: a Package for Flexible Visualization of rearrangement 3D Genome and Multi-omics**
