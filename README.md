[![build status](https://gitlab.quartzbio.com/code/plinker_pkg/badges/master/build.svg)](https://gitlab.quartzbio.com/code/plinker_pkg/commits/master)
[![coverage report](https://gitlab.quartzbio.com/code/plinker_pkg/badges/master/coverage.svg)](https://gitlab.quartzbio.com/code/plinker_pkg/commits/master)

plinker is a R package to interface PLINK over PLINK BED files


## Quick Install

To install the package:

* make install: install the package in the default R library
* sudo make install: if you need root access


# Philosophy

The objective of this package is to use the powerful and fast algorithms implemented
in PLINK without caring about input and ouput file formats, and without relying on PLINK
order of operations.

# Features
  - _plinker_bed_ object:
    * instantaneous loading, .bim file only loaded once if needed
    * is a view of the actual dataset, than can be recursively subsetted
    * a print() method displays relevant information: cardinality, views, annotations...
    
  - subsetting
    * The most powerful feature is to be able to only use subsets of the original dataset without creating new files.
    * you can subset samples or SNPs, by indices or by IDs.
    * you can subset recursively: a subset is an actual first-order object
    * the order of the subsetting is kept: you can reorder your dataset 
      (and the subsequent PLINK algorithms outputs) while keeping the origibnal dataset
    
  - alleles lexigraphic order: the problem with MAF based allele order is that this order may be inversed in a subset.
    * all plinker wrapper PLINK algorithms can use the lexicogaphic order

  - genotypes:
    * seamless and random access to genotypes (numeric) thanks to BEDMatrix R package (N.B: do not need to load .bim)
    * access to genotype strings (using actual alleles)
     
  - custom annotations: You can attach your custom SNP and/or sample annotations to the object, and define new custom IDs.
    If defined, these IDs will appear in the PLINK results.

  - sample IDs:
    * FAM-like annotations are not always well suited: the sample IDs are split between family ID (FID) and 
      internal ID (IID). Quite often, a dummy value is assigned to FID (e.g. "0").
     
    

  - pure-R implementations: for most algorithms, a pure R implementation with same input/output is provided, that allows:
    * to understand what PLINK actually computes
    * to easily tweak it if you need something slighlty different
    * to validate both implementations
     
  - fully tested (coverage ~ 100%)
