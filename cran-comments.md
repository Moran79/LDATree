## Resubmission (08/23/2023)

Sorry for the inconvenience!

This is a resubmission. In this version I have:

* increment the version number from 0.1.0 to 0.1.1

* Fixed two examples in functions due to an MKL error from Additional issues. It detects an out-of-bounds error in one of the plots in the example.

However, the package fails a check on CRAN, and I am not able to understand why. In the error message, it says:

> Pandoc is required to build R Markdown vignettes but not available. Please make sure it is installed.

I do not think I should install "Pandoc" on CRAN's machine. And it passes other flavors except for this r-release-macos-x86_64. Could you please help me with the next step?

Thank you so much!

## Resubmission (08/21/2023)

Thank you for pointing out those problems. Many of these I may never notice.

This is a resubmission. In this version I have:

* Changed all `cat` to verbose style.

* Changed the `<<-` to `<-`

* Successfully remove the archivedFuns.R from the bundle. Used the wrong syntax in the .Rbuildignore file the first time.

About the reference: Sorry, since this method has not been published, I do not have any direct reference to add. Maybe I can do so in future submissions.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
