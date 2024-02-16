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

you can used command line or used api directly.
trackc提供了两种使用模式，1. command line， 2. api

cli可以通过配置文件的方式使用，安装trackc之后，用户可以直接使用trackc命令来调用command line工具
除画图功能之外，trackc集成了常见的Hi-C数据处理命令

可以通过trackc --help来查看

目前集成的功能模块见下
* trackc cli
* trackc gtf2bed
* trackc tadscore
* trackc fa2GC
* trackc compartment
* trackc insulationscore


`Hi-C data analysis guides </analysis_guide/index>`总结了常见Hi-C数据分析教程，可以作为trackc画图的输入文件。



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
