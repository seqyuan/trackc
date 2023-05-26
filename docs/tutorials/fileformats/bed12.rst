=====
gene
=====

The trackc's gene input formats is `BED12 <https://bedtools.readthedocs.io/en/latest/content/general-usage.html#genome-file-format/>`_.

Bed12 file description can be found from the link below: `https://bedtools.readthedocs.io/en/latest/content/general-usage.html#genome-file-format <https://bedtools.readthedocs.io/en/latest/content/general-usage.html#genome-file-format/>`_.

How to generate gene track bed12
=================================

gtf2bed.pl
----------
You can use `gtf2bed` convert gtf format to bed12 format.

`gtf2bed` is a perl script, can be get from the link below: `https://github.com/ExpressionAnalysis/ea-utils/blob/master/clipper/gtf2bed <https://raw.githubusercontent.com/ExpressionAnalysis/ea-utils/master/clipper/gtf2bed>`_.

.. code-block:: shell

    perl gtf2bed GRCh38.84.gtf >GRCh38.84.gtf.bed12

Below table is the output of `gtf2bed`:

.. csv-table:: bed12
   :file: GRCh38.84.gtf.bed12
   :header-rows: 0
   :delim: tab

please note that, the bed12 from gtf2bed.pl is based on transcript id, each row is a transcript, the `column-4` is `transcript id` not gene name


gtf2bed4trackc
--------------
If you have installed `trackc`, you can conver GTF to bed12 using `gtf2bed4trackc` command.
the `column-4` of outfile is `gene name`

.. code-block:: shell

    gtf2bed4trackc -g GRCh38.84.gtf -o GRCh38.84.bed12

.. csv-table:: bed12
   :file: GRCh38.84.bed12
   :header-rows: 0
   :delim: tab
