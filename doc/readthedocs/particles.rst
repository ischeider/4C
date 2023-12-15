Discrete Element Method
=======================

The ``PROBLEM TYP`` particle, viz

::

    -------PROBLEM TYP
    PROBLEMTYP  Particle


enables the simulation of particle problems based on the Discrete Element Method (DEM) via the section ``PARTICLE DYNAMIC``::

    -------PARTICLE DYNAMIC
    INTERACTION     DEM

Stiffness
----------

The algorithm for the stiffness calculation between two particles is defined by the parameter :ref:`NORMALCONTACTLAW<particle_dynamic_dem_normalcontactlaw>`,
given in the section

::

    --------PARTICLE DYNAMIC/DEM
    NORMALCONTACTLAW     <algorithm>

The three most prominent values for ``algorithm`` are
a linear spring model (``NormalLinearSpring``), the Hertzian contact model (``NormalHertz``), and a nonlinear force-displacement curve (``NormalForceDistanceCurve``).
These algorithms may also be extended with damping, which is explained further down.

The force between two particles having a relative normal velocity :math:`v_\text{rel,n}` and a distance 

.. math::
    d = \lVert \mathbf{x}_i - \mathbf{x}_j \rVert_2 - (r_i + r_j) 

(which becomes negative for overlap of interacting particles) is

.. list-table::
    :header-rows: 1

    * - algorithm
      - force displacement relation
    * - NormalLinearSpring
      - :math:`f_n = k_n d`
    * - NormalHertz
      - :math:`f_n = k_n d \sqrt{-d}`
    * - NormalForceDistanceCurve
      - :math:`f_n = f(d)`

If the linear spring or the Hertzian model are used, the normal stiffness :math:`k_n` can either be specified with the input parameter ``NORMAL_STIFF``
or via an estimation based on some (assumed) maximal particle velocity :math:`v_\max` (``MAX_VELOCITY``), density :math:`\rho`,
maximal radius of the involved particles :math:`r_\max` (``MAX_RADIUS``), and a maximum allowed relative (to :math:`r_\max`) penetration :math:`c` (``REL_PENETRATION``).
Then, the stiffness is calculated by

.. math::

    k_n = \frac{2}{3} r_{\max} \Pi \rho \frac{v_{\max}^2}{c^2}.

This calculation assumes a transfer of the complete kinetic energy of two max mass particles approaching centrally into the potential energy of the 
spring acting as repulsive force. As the stiffness :math:`k_n` can be seen as a penalty parameter which hinders the particles to overlap 
it is crucial to have it in a reasonable size, however the exact number is not relevant. Counterproductive large numbers for :math:`k_n` lead to small time step size requirement whereas too small numbers lead to unphysically large overlap of the particles during contact. Note: It is only relevant to keep it constant during
the interaction process of two particles to preserve energy during the impact/restauration process.  

Either the parameter ``NORMAL_STIFF`` or the three parameters ``MAX_VELOCITY, MAX_RADIUS, REL_PENETRATION`` are needed.

In the third option using a force-displacement curve, the relation :math:`f(d)` is given by the function with the number entered as ``NORMAL_CONTACT_FORCE_DIST_CURVE``.

If the force-displacement curve algorithm is used and :math:`c` has a positive value in the input file,
the stiffness in the curve is multiplied by :math:`r_\text{eff}^3`.

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

For the neighborhood search (which strongly influences the performance of the particle problem) two input arguments need to be given
in the binning strategy section

::

    ---------------------BINNING STRATEGY
    DOMAINBOUNDINGBOX  <xmin> <ymin> <zmin> <xmax> <ymax> <zmax>
    BIN_SIZE_LOWER_BOUND <binsize>

The six values in the first row define an axis-aligned cuboid which bounds the region in which particles are tracked. If a particle leaves this region,
it is removed from the simulation. The cuboid is divided into smaller cuboids which have in each direction maximal extension of
``binsize``. The ``binsize`` should be chosen such that it is larger than twice the maximal travel distance of a particle within a time step added by the
interaction range of particles (e.g. their radii). In doing so, no collision event is missed. An error is thrown if this criteria is violated.
Either the time step size can be reduced or the ``binsize`` can be enlarged to resolve the issue.

Instead of specifying ``BIN_SIZE_LOWER_BOUND`` it is also possible to give the number of bins per spatial direction via

::

    BIN_PER_DIR <num_x_dir> <num_y_dir> <num_z_dir>

such that the big cuboid is subdivided. 

Note: The bins are a used for parallel distribution of the problem when  a multi-core simulation is performed.
 
If a periodic box of particles is necessary, the following line can be added in the binning strategy section

::

    PERIODICONOFF    1 1 1

The three entries' default value is ``0`` meaning no peridicity in the corresponding axis direction (:math:`x`, :math:`y`, :math:`z`).
A ``1`` activates the periodicity in the respective spatial direction.

Note that the the extension of the :ref:`DOMAINBOUNDINGBOX<binningstrategy_domainboundingbox>`
must be at least three times of :ref:`BIN_SIZE_LOWER_BOUND<binningstrategy_bin_size_lower_bound>` in a periodic direction.


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
