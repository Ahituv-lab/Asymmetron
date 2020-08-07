.. _orientation.rst:
  
===========
orientation
===========

|

**This function performs strand assignment of an unassigned feature based on another overlapping feature, thereby enabling the strand asymmetry analysis of the first.**

|

.. list-table:: **Required inputs:**
   :header-rows: 1

   * - Argument
     - Explanation

   * - motif_no_annotation
     - Un-annotated file to which strand assignment will be added.
   * - motif_annotation
     - Strand-annotated file used to extract strand annotation from.


|


**Output:**

* Un-annotated BED formatted motif with annotation based on motif_annotation file.

