MPN v0.4.0
============

Changes:

* Added a 'tol' argument to mpn() and apc(). This value will be passed to
stats::uniroot() and stats::optimize(). The default value is lower than
previously used.
* Updated URL for the associated shiny app.
* Minor updates to documentation.


MPN v0.3.0
============

Changes:

* Added the apc() function to estimate Aerobic Plate Counts (APC).
* Added a bias-adjusted point estimate to the mpn() output.


MPN v0.2.0
============

Changes:

* Added a 'CI_method' argument with a likelihood ratio confidence interval
option to mpn().
* Minor adjustments to uniroot() argument values.
