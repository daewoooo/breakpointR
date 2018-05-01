![breakpointR](breakpointR_logo.png)
====================================

An R-package for breakpoint detection in single cell Strand-seq data.

Collaborators: David Porubsky, Ashley D Sanders, Aaron Taudt

Installation
------------

### Stable bioconductor version (not available yet)
Under development.

### Development version from Github
To install the development version from Github, follow the steps given below. The installation has only been tested on Ubuntu so far, if you need to install on Windows or Mac additional steps might be necessary (e.g. installation of Rtools from https://cran.r-project.org/bin/windows/Rtools/)

1. Install a recent version of R (>=3.4.0) from https://www.r-project.org/
2. Optional: For ease of use, install Rstudio from https://www.rstudio.com/
3. Open R and install all dependencies. Please ensure that you have writing permissions to install packages. Execute the following lines one by one:

   install.packages("devtools")  
	 source("http://bioconductor.org/biocLite.R")  
	 biocLite("GenomicRanges")  
	 biocLite("GenomicAlignments")  
	 library(devtools)
	 install_github("daewoooo/strandseqExampleData") 
	 install_github("daewoooo/BreakPointR")  
	 #### Or alternatively if the above line doesn't work:
	 install_git("git://github.com/daewoooo/strandseqExampleData.git", branch = "master")
	 install_git("git://github.com/daewoooo/BreakPointR.git", branch = "master")

How to use BreakPointR
----------------------

Please refert to the [vignette](https://github.com/daewoooo/BreakPointR/blob/master/vignettes/breakpointR.pdf) for tutorials on all breakpointR features.

Report Errors
-------------

If you encounter errors of any kind, please report an [issue here](https://github.com/daewoooo/BreakPointR/issues/new).
