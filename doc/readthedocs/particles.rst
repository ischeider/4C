Discrete Element Method
=======================

The Discrete Element Method (DEM) is defined by the section ``PARTICLE DYNAMIC``::

    -------PARTICLE DYNAMIC
    INTERACTION     DEM

Stiffness
----------

The algorithm for the stiffness calculation between two particles is defined by the parameter :ref:`NORMALCONTACTLAW<particle_dynamic_dem_normalcontactlaw>`,
given in the section

::

    --------PARTICLE DYNAMIC/DEM
    NORMALCONTACTLAW     <algorithm>

The two most prominent values for ``algorithm`` are
a linear spring model (``NormalLinearSpring``), the Hertzian contact model (``NormalHertz``) and a nonlinear force-displacement curve (``NormalForceDistanceCurve``).
These algorithms may also be extended with damping, which is explained further down.

The force between two particles having a relative normal velocity :math:`v_\text{rel,n}` and a distance :math:`d` is

.. list-table::
    :header-rows: 1

    * - algorithm
      - force displacement relation
    * - NormalLinearSpring
      - :math:`f_n = K_n d`
    * - NormalForceDistanceCurve
      - :math:`f_n = f(d)`
    * - NormalHertz
      - :math:`f_n = K_n d \sqrt{-d}`

If the linear spring or the Hertzian model are used, :math:`K_n` is given by the parameter ``NORMAL_STIFF``.
For the force-displacement curve, the relation :math:`f(d)` is given by the function with the number entered as ``NORMAL_CONTACT_FORCE_DIST_CURVE``.
Another available option is the calculation of the stiffness
from some (assumed) maximum velocity :math:`v_\max` and a maximum allowed penetration, :math:`c`.
Then, the stiffness is calculated by

.. math::

    K_n = \frac{m_P \left( \frac{v_\max}{c}\right)^2}{A_p}

with

.. math::

    \begin{split}
    m_P &= \frac{4}{3} \rho \pi r_\text{eff}^3 \\
    A_p &= 2 r_\text{eff}^2 \\
    r_\text{eff} &= 0.5 (r_i + r_j)
    \end {split}

That is, the stiffness is linear in the radius of the particle.

If the force-displacement curve algorithm is used and :math:`c` has a positive value in the input file,
the stiffness in the curve is multiplied by :math:`r_\text{eff}^3`.

The values for :math:`v_\max, c` correspond to the damping parameters ``MAX_VELOCITY, REL_PENETRATION``.

Thus, either the parameter ``NORMAL_STIFF`` or the two parameters ``MAX_VELOCITY, REL_PENETRATION`` are needed.

.. note::

    *Personal notes:

    #. Actually, I don't know why the area in the stiffness calculation is calculated this way.
    #. If the stiffness is given by a force-displacement curve, the force is rather a force-per-volume value.
    #. Previously, the calculation of the stiffness depending on an artificial mass did not make much sense, if the mass changes with time
       (A constant parameter ``MAX_RADIUS`` was used, and thus the stiffness was constant).
       This has been modified, and ``MAX_RADIUS`` is not needed anymore.

Damping
--------

Damping can be applied in two ways: First, there is a viscous damping algorithm, which acts as an artificial damping matrix to damp the whole system.
This can be switched on by a value greater 0 for the parameter ``VISCOUS_DAMPING``::

    -------PARTICLE DYNAMIC
    VISCOUS_DAMPING     0.1

The second option is to modify the stiffness of the particle interaction.
Here, damping is included in the following algorithms:

- NormalLinearSpringDamp
- NormalForceDistanceCurveDamp
- NormalLeeHerrmann
- NormalKuwabaraKono
- NormalTsuji

a subset of the following parameters

::

    ----------PARTICLE DYNAMIC/DEM
    COEFF_RESTITUTION
    DAMP_REG_FAC
    NORMAL_DAMP

may be given in the input file, depending on the algorithm. In general the force-distance curve is given by the following equations

.. list-table::
    :header-rows: 1

    * - algorithm
      - force displacement relation
    * - NormalLinearSpringDamp
      - :math:`f_n = K_n d - d_n v_\text{rel,n} f_\text{reg}`
    * - NormalForceDistanceCurveDamp
      - :math:`f_n = f(d) - d_n v_\text{rel,n} f_\text{reg}`
    * - NormalLeeHerrmann
      - :math:`f_n = K_n d \sqrt{-d} - m_\text{eff} d_n v_\text{rel,n}`
    * - NormalKuwabaraKono
      - :math:`f_n = K_n d - d_n v_\text{rel,n} \sqrt{-d};`
    * - NormalTsuji
      - :math:`f_n = K_n d - d_n v_\text{rel,n} \sqrt{\sqrt{-d}};`

In the three latter models, the damping factor, :math:`d_n` is given by the input parameter ``NORMAL_DAMP``.

Some damping models depend on the effective mass of the two interacting particles, :math:`m_\text{eff}`.
This is calculated based on the current radius and density by

.. math::

    m_\text{eff} = \frac{m_1 m_2}{m_1 + m_2}


For the linear spring and the force-distance curve models, the damping factor :math:`d_n` is calculated depending on the value :math:`e` (corresponding to input parameter ``COEFF_RESTITUTION``) by

.. math::

    d_n = \left\{
    \begin{array}{cl} 2 \sqrt{K_n m_\text{eff}} & e=0 \text{ or unset} \\
        2 \sqrt{K_n m_\text{eff}} \frac{|\ln{e}|}{\sqrt{\ln{e}^2 + \pi^2}} & e>0
    \end{array}\right.

The regularization factor :math:`f_\text{reg}` is either 1, if the value :math:`c_{r_\text{fac}}` corresponding to the input parameter ``DAMP_REG_FAC`` is unset (or 0),
or it is calculated depending on the distance and the current radii of the two particles:

.. math::

    f_\text{reg} = \min\left(1,\, \frac{d}{c_{r_\text{fac}} \min(r_1,r_2)}\right)



Summary and advice
-------------------

Summarizing the above, the simplest way for an undamped particle-to-particle relation can be produced by

::

    --------PARTICLE DYNAMIC/DEM
    NORMALCONTACTLAW     NormalLinearSpring
    NORMAL_STIFF         <value>

A bit more complicated is the use of a maximum velocity and penetration instead of the stiffness,
which then makes the stiffness depending on the radius (and thus the mass) of the particles.

If the particles change their size during the simulation, one should not prescribe the stiffness with a constant value,
but depending on the velocity and the maximum allowed penetration.

If damping is used, the least number of parameters is used
if :math:`c_{f_\text{fac}}` and :math:`e` (parameters ``DAMP_REG_FAC`` and ``COEFF_RESTITUTION``) remain unset.
For this case, no further parameter is needed, since the damping factor :math:`d_n` is then calculated by the stiffness and the mass of the particles.

If  :math:`e>0`, it is included in the calculation of the damping factor :math:`d_n`.

If additionally (or even alternatively) the factor :math:`c_{f_\text{fac}} > 0`,
the regularization factor depends linearly on the ratio of distance to radius.

Bounding Box
-------------

If one would like to have a periodic box of particles (like atoms),
one has to describe the following in the binning strategy::

    ---------------------BINNING STRATEGY
    PERIODICONOFF    1 1 1
    DOMAINBOUNDINGBOX  <xmin> <ymin> <zmin> <xmax> <ymax> <zmax>
    BIN_SIZE_LOWER_BOUND <binsize>

Note that the the extension of the :ref:`DOMAINBOUNDINGBOX<binningstrategy_domainboundingbox>`
must be at least three times of :ref:`BIN_SIZE_LOWER_BOUND<binningstrategy_bin_size_lower_bound>` in every direction.


Output
-------

Direct runtime output is available for all particle related result values.

If output of the bounding box is intended, add the keyword ``WRITEBINS`` to the section ``BINNING STRATEGY``::

    ---------------------BINNING STRATEGY
    WRITEBINS    rows|cols

``rows`` produces HEX27 elements for all bins that contain the owned particles,
while ``cols`` also includes the ghosted particles, that is, the HEX elements will overlap
(which is usually not what one might want).

If the hex elements shall be visualized within Paraview, one has to run the post_processor tool, so that the hex elements are transferred into a format that can be read by Paraview.

::

    post_processor --file=<simulation_resultbase>

is sufficient to get an ensight file ``<simulation_resultbase>_bins.case``
