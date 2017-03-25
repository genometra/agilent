
agilent
================================================================================

[![Build Status](https://travis-ci.org/genometra/agilent.svg?branch=master)](https://travis-ci.org/genometra/agilent)

An R library to normalize Agilent one color microarray data.

Includes functions for:

- parsing Agilent headers
- reading expression, copy number or arrayCGH data
- reading data into an expression set
- reading Agilent or GPR scanned arrays
- background correction
- normalization
- average duplicated spots



Install
--------------------------------------------------------------------------------

The latest version of the library can be installed form the R session doing:

    install.packages ("devtools")
    library (devtools)
    setRepositories (ind=1:2)
    install_github ("genometra/agilent/pkg")
