Overview
--------

archetypal is a package for performing Archetypal Analysis (AA) by using
a properly modified version of PCHA algorithm.

Basic functions are:

-   `archetypal()` do AA
-   `find_outmost_projected_convexhull_points` Projected CH initial
    solution.
-   `find_outmost_convexhull_points` CH initial solution.
-   `find_outmost_partitioned_convexhull_points()` Partitioned CH
    initial solution.
-   `find_furthestsum_points()` Furthest Sum initial solution.
-   `find_outmost_points()` Outmost initial solution.
-   `find_optimal_kappas()` search for the optimal number of archetypes
-   `find_pcha_optimal_parameters()` search for the optimal updating
    parameters of PCHA algorithm 
-   `check_Bmatrix()` check B matrix after run of AA.
-   `study_AAconvergence()` study the convergence of PCHA algorithm
-   `find_closer_points()` find the closer to archetypes data points

Install the archetypal package and then read
`vignette("archetypal", package = "archetypal")`.

Installation
------------

``` r
# Install with dependencies:
install.packages("archetypal",dependencies=TRUE)
```

Usage
-----

``` r
library(archetypal)

data("wd2")
df = wd2
aa = archetypal(df = df, kappas = 3,verbose = FALSE, rseed = 9102)

# Time for computing Projected Convex Hull was 0.01 secs 
# Next projected convex hull initial solution will be used... 
#           x        y
# 34 5.687791 3.481611
# 62 1.961799 2.793497
# 5  5.123878 2.745874
# 
# archs=aa$BY
# archs
# x        y
# [1,] 5.430757 3.146258
# [2,] 2.043435 2.710947
# [3,] 3.128401 4.781751
# aa[c("SSE","varexpl","iterations","time" )]
# $SSE
# [1] 1.717538
# 
# $varexpl
# [1] 0.9993186
# 
# $iterations
# [1] 63
# 
# $time
# [1] 8.1
# cbind(names(aa))
# [,1]             
# [1,] "BY"             
# [2,] "A"              
# [3,] "B"              
# [4,] "SSE"            
# [5,] "varexpl"        
# [6,] "initialsolution"
# [7,] "freqstable"     
# [8,] "iterations"     
# [9,] "time"           
# [10,] "converges"      
# [11,] "nAup"           
# [12,] "nAdown"         
# [13,] "nBup"           
# [14,] "nBdown"         
# [15,] "run_results"   
```

Contact
-------

Please send comments, suggestions or bug reports to
<dchristop@econ.uoa.gr> or <dem.christop@gmail.com>