This version is submitted because the accompanying paper is about to be released
by the Journal of Statistical Software. They have asked that we update the DOI
and citation details in CITATION. This version also fixes minor bugs.

## R CMD check results

0 errors | 0 warnings | 1 note

Note:

```
Found the following (possibly) invalid DOIs:
     DOI: 10.18637/jss.v115.i02
       From: DESCRIPTION
             inst/CITATION
             man/sdmTMB.Rd
             man/sdmTMB.Rd
       Status: 404
       Message: Not Found
```

This is expected. The DOI will not be active until the journal submits the DOI,
which will happen after this package update on CRAN.

## Test environments

* Local M2 macOS install, R 4.5.2
* Intel macOS (on github-actions), R-release
* Ubuntu 24.04.3 (on github-actions), R-release
* Windows (on github-actions), R-release
* Windows (winbuilder), R-devel
* Windows (winbuilder), R-release
* Windows (winbuilder), R-oldrelease
