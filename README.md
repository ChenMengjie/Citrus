### Citrus

Citrus is a toolkit for single cell sequencing analysis under development. 
Citrus is written in c++ with easy-to-use R wrappers. 
To use Citrus, the installation of **Rcpp** and **RcppArmadillo** is required. 

### Available methods

Now a normalized method and a clustering method are available in the distributed version. 

* Method 1: **scPLS**: a normalization method to remove unwanted variation using both control and target genes.

  The description of methodology and algorithm can be found in the manuscript:

  Chen M, Zhou X: **Normalization of single cell RNA sequencing data using both control and target genes**. 2016

  which can be downloaded at: http://biorxiv.org/content/early/2016/03/21/045070.

  A **vignette** on how to use scPLS is distributed together with the package.

* Method 2: a clustering method which can account for confounding factors.

  The description of method and a tutorial of software on clustering will be available shortly.


### Installation

**Citrus** relies on the following R packages: **Rcpp**, **RcppArmadillo**, **genlasso**, **gplots**. All packagess are hosted on CRAN. 
  ```R
  install.packages ("Rcpp")
  install.packages ("RcppArmadillo")
  install.packages ("genlasso")
  install.packages ("gplots")
  ```

**Citrus** can be installed from github directly as follows:

  ```R
  install.packages ("devtools")
  library(devtools)
  install_github("ChenMengjie/Citrus")
  ```
  
### Why do we call it citrus? 
Citrus is not a short name for any of our methods. To honor Bowtie/Cufflinks and Salmon/Sailfish series, we are creating a fruit collection. 

### Authors

**Mengjie Chen** (U Chicago) 
**Xiang Zhou** (Umich)

Bug report, comments or questions please send to mengjiechen@uchicago.edu.
