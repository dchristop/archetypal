# archetypal 1.1.1 (2020-10-09)

* fix a bug for OS r-patched-solaris-x86


# archetypal 1.1.0 (2020-01-27) - major update

## New functions

* `find_pcha_optimal_parameters()` finds the optimal parameters to be used for PCHA algorithm. 
* `find_closer_points()` finds the data points that are closer to the archetypes
during all iterations of algorithm PCHA.
* `study_AAconvergence()` studies the convergence of Archetypal Analysis when
using algorithm PCHA.
* `grouped_resample()` performs simple or Dirichlet resampling.
* `dirichlet_sample()` performs Dirichlet sampling.

## New data sets

* `Absolute Temperature` the Global Absolute Temperature data set for Northern Hemisphere 1969-
2013.

* `gallupGPS6` the Gallup Global Preferences Study processed data set of six variables.

* `wd25` a 2D data set created by 5 points for demonstration purposes.

## Minor changes

* `find_outmost_projected_convexhull_points` changes its `n` argument to `npr` and has new arguments `rseed`, `doparallel`, `nworkers` and 
`uniquerows`.
* `check_Bmatrix` changes its `print.details` argument to `verbose`.
* `find_furthestsum_points` changes its `nworkers = 10` argument to `nworkers = NULL` and has new argument  `doparallel`.
* `align_archetypes_from_list` has a new argument `verbose`.

