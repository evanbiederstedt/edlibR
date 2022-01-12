[![<evanbiederstedt>](https://circleci.com/gh/evanbiederstedt/edlibR.svg?style=svg)](https://app.circleci.com/pipelines/github/evanbiederstedt/edlibR)
[![CRAN status](https://www.r-pkg.org/badges/version/edlibR)](https://cran.r-project.org/package=edlibR)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/edlibR)](https://cran.r-project.org/package=edlibR)

 
# edlibR

- [Introduction](#edlibr-r-integration-for-edlib)
- [Main Features](#main-features)
- [Usage](#usage)
- [Background](#background)
  * [Alignment Methods](#alignment-methods)
  * [Algorithmic Design](#algorithmic-design)
- [Installation](#installation)
- [References](#references)

## edlibR: R integration for edlib

The R package `edlibR` provides bindings to [edlib](https://github.com/Martinsos/edlib), a performant C/C++ library for calculating the exact pairwise sequence alignment using the [edit distance](https://en.wikipedia.org/wiki/Edit_distance) ([Levenshtein distance](https://en.wikipedia.org/wiki/Levenshtein_distance)). The package is modeled after the API of the [Python package edlib](https://pypi.org/project/edlib/).

The original Python bindings to the C/C++ library using Cython can be found [here](https://github.com/Martinsos/edlib/tree/master/bindings/python). Much of the background information below is directly from the [original README for edlib](https://github.com/Martinsos/edlib/blob/master/README.md).

## Main Features

* It calculates the **edit distance (Levenshtein distance)**.
* It can find the **optimal alignment path** (instructions how to transform the first sequence into the second sequence).
* It can find just the **start and/or end locations of the alignment path** - can be useful when speed is more important than having the exact alignment path.
* Supports **multiple alignment methods**: global(**NW**), prefix(**SHW**) and infix(**HW**), each of them useful for different scenarios.
* You can **extend character equality definition**, enabling you to e.g. have wildcard characters, to have case insensitive alignment or to work with degenerate nucleotides.
* It can easily handle small or **very large sequences**, even when finding the alignment path, while consuming very little memory.
* **Super fast** thanks to Myers's bit-vector algorithm, as described in ["A Fast Bit-Vector Algorithm for Approximate String Matching Based on Dynamic Programming"](http://www.gersteinlab.org/courses/452/09-spring/pdf/Myers.pdf).

## Usage

For details and examples of the edlibR functionality, see the **vignettes:**
* [HTML version](https://htmlpreview.github.io/?https://raw.githubusercontent.com/evanbiederstedt/edlibR/blob/main/doc/edlibr_walkthrough.html)
* [Markdown version](https://github.com/evanbiederstedt/edlibR/blob/main/vignettes/edlibR_walkthrough.md)

There are three functions within edlibR:

**align()**

```
align(query, target, [mode], [task], [k], [cigarFormat], [additionalEqualities])
```

**getNiceAlignment()**

```
getNiceAlignment(alignResult, query, target, [gapSymbol])
```

**nice_print()**

```
nice_print(niceAlignment)
```

## Background

### Alignment Methods

Edlib supports three alignment methods:
* **global (NW)** - This is the standard method, when we say "edit distance" this is the method that is assumed. (Note that 'NW' stands for 'Needleman-Wunsch'.)
  It tells us the smallest number of operations needed to transform the first sequence into the second sequence.
  *This method is appropriate when you want to find out how similar the first sequence is to the second sequence.*
* **prefix (SHW)** - Similar to the global method, but with a small twist - the gap at the query end is not penalized. What this means is that deleting elements from the end of the second sequence is "free"! (Note that 'HW' stands for 'Hybrid Wunsch'.)
  For example, if we had query as `AACT` and target as `AACTGGC`, the edit distance would be 0, because removing `GGC` from the end of the second sequence is "free" and does not count into the total edit distance.
  *This method is appropriate when you want to find out how well the first sequence fits at the beginning of the second sequence.*
* **infix (HW)**: Similar as prefix method, but with one more twist - gaps at query end **and start** are not penalized. What this means is that deleting elements from both the start and the end of the second sequence is "free"! (Note that 'SHW' stands for 'Semi-Hybrid Wunsch'.)
  For example, if we had query as `ACT` and target as `CGACTGAC`, the edit distance would be 0, because removing `CG` from the start and `GAC` from the end of the second sequence is "free" and does not count into the total edit distance.
  *This method is appropriate when you want to find out how well the first sequence fits at any part of the second sequence.* For example, if your second sequence was a long text and your first sequence was a sentence from that text, but slightly scrambled, you could use this method to discover how scrambled it is and where it fits in that text.
  *In bioinformatics, this method is appropriate for aligning a read to a sequence.*


### Algorithmic Design

Edlib is based on Myers's bit-vector algorithm[[1]](#1) and extends from it.

It calculates a dynamic programming matrix of dimensions `Q x T`, where `Q` is the length of the first sequence (query), and `T` is the length of the second sequence (target). It uses Ukkonen's banded algorithm[[2]](#2) to reduce the space of search, and there is also parallelization from Myers's algorithm, however time complexity is still quadratic.

Edlib uses Hirschberg's algorithm[[3]](#3) to find the alignment path, therefore space complexity is linear.

More details can be found within the edlib publication[[4]](#4).


## Installation
 
 
To install the stable version from [CRAN](https://cran.r-project.org/package=edlibR), use:

```r
install.packages('edlibR')
```


To install the latest version, use:

```r
install.packages('devtools')
devtools::install_github('evanbiederstedt/edlibR')
```


## References

### R package citation

The R package can be cited as follows:

```
Martin Šošić and Evan Biederstedt (2022). edlibR: R wrapper to edlib,
the C/C++ library for fast exact sequence alignment using edit
(Levenshtein) distance. R package version 1.0.0.
https://github.com/evanbiederstedt/edlibR
```

### edlib publication

Please also cite the original edlib publication, which can be found [here](https://academic.oup.com/bioinformatics/article/33/9/1394/2964763) and should be cited as follows:

```
Martin Šošić, Mile Šikić (2017). Edlib: a C/C++ library for fast, exact sequence alignment using edit distance. 
Bioinformatics, Volume 33, Issue 9, 1 May 2017, Pages 1394–1395, 
https://doi.org/10.1093/bioinformatics/btw753
```

### README references

<a id="1">[1]</a> 
Myers G. (1999) A fast bit-vector algorithm for approximate string matching based on dynamic programming. J. ACM, 46, 395–415. https://doi.org/10.1145/316542.316550

<a id="2">[2]</a>
Ukkonen E. (1985) Algorithms for approximate string matching. Inform. Control, 64, 100–118. https://doi.org/10.1016/S0019-9958(85)80046-2

<a id="3">[3]</a>
Hirschberg D.S. (1975) A linear space algorithm for computing maximal common subsequences. Commun. ACM, 18, 341–343. https://doi.org/10.1145/360825.360861

<a id="4">[4]</a>
Martin Šošić, Mile Šikić. Edlib: a C/C++ library for fast, exact sequence alignment using edit distance, Bioinformatics, Volume 33, Issue 9, 1 May 2017, Pages 1394–1395, 
https://doi.org/10.1093/bioinformatics/btw753
