			R-Forge SVN README

This file explains the repository structure of your project. A more
detailed guide to R-Forge is available by 
Theußl and Zeileis (2010) [1] and the R-Forge Administration and 
Development Team (2009) [2].

1. Introduction
-----------------------------------------------------------------------
R is free software distributed under a GNU-style copyleft. R-Forge is
a central platform for the development of R packages, R-related 
software and further projects. Among many other web-based features it 
provides facilities for collaborative source code management via 
Subversion (SVN) [3].

2. The directory you're in
-----------------------------------------------------------------------
This is the repository of your project. It contains two important
pre-defined directories namely 'pkg' and 'www'. These directories must 
not be deleted otherwise R-Forge's core functionality will not be 
available (i.e., daily checking and building of your package or the 
project websites).
'pkg' and 'www' are standardized and therefore are going to be
described in this README. The rest of your repository can be used as
you like.

3. 'pkg' directory
-----------------------------------------------------------------------
To make use of the package building and checking feature the package 
source code has to be put into the 'pkg' directory of your repository 
(i.e., 'pkg/DESCRIPTION', 'pkg/R', 'pkg/man', etc.) or, alternatively,
a subdirectory of 'pkg'. The latter structure allows for having more 
than one package in a single project, e.g., if a project consists of 
the packages foo and bar then the source code will be located in 
'pkg/foo' and 'pkg/bar', respectively.

R-Forge automatically examines the 'pkg' directory of every repository 
and builds the package sources as well as the package binaries on a
daily basis for Mac OS X and Windows (if applicable). The package builds
are provided in the 'R Packages' tab for download or can be installed
directly in R from a CRAN-style repository using 
'install.packages("foo", repos="http://R-Forge.R-project.org")'. 
Furthermore, in the 'R Packages' tab developers can examine logs 
generated on different platforms by the build and check process.

4. 'www' directory
-----------------------------------------------------------------------
Developers may present their project on a subdomain of R-Forge, e.g.,
'http://foo.R-Forge.R-project.org', or via a link to an external
website.

This directory contains the project homepage which gets updated hourly
on R-Forge, so please take into consideration that it will not be 
available right after you commit your changes or additions. 

5. Help
-----------------------------------------------------------------------
If you need help don't hesitate to submit a support request at 
https://r-forge.r-project.org/tracker/?func=add&group_id=34&atid=194, 
search the forum 
https://r-forge.r-project.org/forum/forum.php?forum_id=78&group_id=34,
or contact us at R-Forge@R-project.org.

6. References
-----------------------------------------------------------------------

[1] Stefan Theußl and Achim Zeileis. Collaborative software development 
using R-Forge. The R Journal, 1(1):9-14, May 2009. URL 
http://journal.r-project.org/2009-1/RJournal_2009-1_Theussl+Zeileis.pdf  

[2] R-Forge Administration and Development Team. RForge User’s Manual, 
2008. URL http://download.R-Forge.R-project.org/R-Forge.pdf

[3] C. M. Pilato, B. Collins-Sussman, and B. W. Fitzpatrick. Version 
Control with Subversion. O’Reilly, 2004. Full book available online at 
http://svnbook.red-bean.com/
