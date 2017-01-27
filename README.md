# IDR(s) Solver

The Induced Dimension Reduction method (IDR(s)) [1] is a short-recurrences Krylov method that
solves the system of linear equation,
                                      Ax = b.
This c++ implementation uses the Matrix Templates Library [MTL4](http://www.simunova.com/home) and is based on [2]. 

References
----------

    .. [1] P. Sonneveld and M. B. van Gijzen
             SIAM J. Sci. Comput. Vol. 31, No. 2, pp. 1035--1062, (2008).
    .. [2] M. B. van Gijzen and P. Sonneveld
             ACM Trans. Math. Software,, Vol. 38, No. 1, pp. 5:1-5:19, (2011).
    .. [3] This file is a translation of the following MATLAB implementation:
            http://ta.twi.tudelft.nl/nw/users/gijzen/idrs.m
