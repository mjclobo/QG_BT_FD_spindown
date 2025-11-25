This is a doubly-reentrant barotropic QG model
===================================

.. note::

   This project is under active development.

Introduction
--------
This is a publically available barotropic quasi-geostrophic model that
can be used for pedagogy or research.

In its present form, the model replicates the results of Nakamura and Plumb (1994)
by using time-varying topography to simulate wave breaking in Earth's atmosphere.
Instead of using the contour advection method of Nakamura and Plumb (1994),
we use finite-difference methods (it took me an embarrassingly long time to
conclude that spectral methods can't capture the sharp jump in the initial conditions).

Any questions or corrections, whether related to the program or related to this documentation, may be
sent to *mattlobo@princeton.edu*.

Contents
--------

.. toctree::

   modeleqns
   exampleruns
   
