# HVCF

HVCF is a pilot project to evaluate pros/cons of storing VCF using HDF5.

Features:
* Variants are indexed by name (using hash-based index stored on disk) and position (using interval-based index stored on disk).
* Samples are indexed by name (using hash-based index stored on disk).
* Efficient genotypes/haplotypes subsetting by position and sample (i.e. row and column) simultaneously.
* Supports compression using GZIP and Blasc LZ4HC.
* Build-in LD (r, r^2) computation using Armadillo linear algebra library.
