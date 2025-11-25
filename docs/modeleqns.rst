=========================================
Model equations
=========================================
Here we provide that governing equations of a barotropic QG model on a doubly-reentrant
:math:`f` plane.

Governing equations
~~~~~~~~~~~~~~~~~~~~~~~~~
The governing equations for the present two-layer QG model are

.. math::

   \begin{equation}
   \partial_{t} q
   + \mathrm{J} (\psi_{v} + \psi_{f}, \, q )
   = 
   - \nu \nabla^4 q
   \end{equation}

where the potential vorticity is
.. math::

    q = f + \nabla^2 \psi - \lambda^{-2} \psi + f \, h / H

and we use a second-order
hyperviscosity with the coefficient :math:`\nu` with dimensions of :math:`\mathrm{L}^{4} \mathrm{T}`.

The time-dependent topographic forcing is defined as

.. math::

   \psi_{f} (x,y,t)
   =
   \psi_{0} \mathrm{sin} ( 2 \pi t / \tau ) \times
   \mathrm{exp} \left [
   - \left {
   \frac{(x-x_{0})^2}{2 \, \delta x^2}
   +  \frac{(y-y_{0})^2}{2 \, \delta y^2}
   \right }
   \right ]

The velocity streamfunction is defined via
.. math::

    (\nabla^2 - \lambda^{-2}) \psi_{v} = q(x, y, t) - f

Note that the PV is advected by the total streamfunction,
:math"`\psi_{f} + \psi_{v}`.

Refer to Nakamura and Plumb (1994) for more details.
Though default parameters reproduce their results, all parameters are easily
modified, as we now show on the next page.



