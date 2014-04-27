# mgarch - R Package for MGARCH Model Simulation and Estimation

*mgarch* is an R package for simulating and estimating multivariate
GARCH models. Currently, only BEKK, GJR and DCC specifications are
implemented.

## Basic Usage

After you install the package, you can load the library:

    > library(mgarch)

Probably, you are interested in the following functions:

1. `mvBEKK.sim`: Simulates an MGARCH-BEKK process
2. `mvBEKK.est`: Estimates an MGARCH-BEKK model
3. `mvBEKK.diag`: Provides diagnostics for an estimated MGARCH-BEKK
   model

You can get help from the embedded documentation by issuing:

    > ?mvBEKK.sim
    > ?mvBEKK.est
    > ?mvBEKK.diag

As an example, you can simulate an MGARCH-BEKK process as follows:

    > sim <- mvBEKK.sim(series.count=2, T=1200)

This will give you two vectors of 2500 observations along with many
other objects. Checkout the help for further details.

Now create a data frame with the given two vectors:

    > eps <- data.frame(s1=sim$eps[[1]], s2=sim$eps[[2]])

The first and last few rows of the resulting data frame looks like as
follows:

    > head(eps)
           s1      s2
    1 -0.7302 -0.9066
    2 -1.0601 -0.9100
    3 -2.0181  0.8561
    4  1.1678  0.9071
    5 -0.1927  0.2441
    6  0.9193 -0.0993

    > tail(eps)
              s1      s2
    2495 -1.0160 -0.9626
    2496  1.5960 -0.5020
    2497 -1.6413  1.1664
    2498  2.0670 -1.9811
    2499  0.6587  0.1838
    2500  1.0310  1.2592

When you need to supply your own data, make sure that it follows the
same convention as above.

Now, we will estimate a model:

    > est <- mvBEKK.est(eps)

The estimation may take longer depending on the number of series and
the length of series. The result contains quite a lot of information:

    > names(est)
    [1] "eps"               "series.length"
    [3] "estimation.time"   "total.time"
    [5] "order"             "estimation"
    [7] "aic"               "asy.se.coef"
    [9] "est.params"        "cor"
    [11] "sd"                "H.estimated"
    [13] "eigenvalues"       "uncond.cov.matrix"
    [15] "residuals"

Check the documentation for details.

You can see a summary by using the `mvBEKK.diag` command:

    > mvBEKK.diag(est)

        Number of estimated series :  2
        Length of estimated series :  1200
        Estimation Time            :  3.971
        Total Time                 :  4.788
        BEKK order                 :  1 1
        Eigenvalues                :  0.9371 0.5384 0.5016 0.3784
        aic                        :  5592
        unconditional cov. matrix  :  3.821 -0.5761 -0.5761 11.3
        var(resid 1 )                :  0.9996
        mean(resid 1 )               :  0.0005716
        var(resid 2 )                :  0.9922
        mean(resid 2 )               :  0.04112
        Estimated parameters       :

        C estimates:
         [,1]   [,2]
    [1,] 1.16 0.2905
    [2,] 0.00 0.9404

        ARCH estimates:
            [,1]     [,2]
    [1,] -0.5225  0.14426
    [2,]  0.3576 -0.09822

        GARCH estimates:
            [,1]     [,2]
    [1,] -0.5358 -0.07182
    [2,] -0.2370  0.90391

        asy.se.coef                :

        C estimates, standard errors:
           [,1]   [,2]
    [1,] 0.1767 0.7308
    [2,] 0.0000 0.2853

        ARCH estimates, standard errors:
            [,1]    [,2]
    [1,] 0.03618 0.03759
    [2,] 0.03326 0.03972

        GARCH estimates, standard errors:
           [,1]    [,2]
    [1,] 0.0519 0.09403
    [2,] 0.1336 0.03051
    Called from: mvBEKK.diag(est)
    Browse[1]>

This command will also produce some helpful plots.

## Copyrights and Licensing

This software is licensed under the GPLv3::

    Copyright (C) 2005-2014 Harald Schmidbauer <harald@hs-stat.com>,
        Angi Roesch <angi@angi-stat.com>,
        Vehbi Sinan Tunalioglu <vst@vsthost.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
