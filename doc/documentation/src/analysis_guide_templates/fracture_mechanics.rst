.. _fracture_mechanics:

Fracture mechanics
==================

4C can evaluate a scalar J-integral along a three-dimensional crack front using an equivalent
domain integral (EDI). Results are evaluated independently for every crack-front node and for
multiple topology-based contour layers. The implementation also retains the configurational-force
vector and reports its projection onto the plane normal to the crack-front tangent as a predicted
propagation direction. It does not grow the crack.

Input
-----

The crack front is identified by an external mesh node-set name or by the one-based ID of a legacy
``DLINE`` node set:

.. code-block:: yaml

   STRUCTURAL DYNAMIC/FRACTURE:
     - crack_id: crackfront
       direction: [1.0, 0.0, 0.0]
       normal: [0.0, 1.0, 0.0]
       contours: 5
       propagation_criterion: maximum_configurational_force
       output:
         quantity: [J, G, propagation_direction]
         output_contours: all
         data_format: ascii

The crack-front nodes must coincide with one connected chain or cycle of solid-element edges.
``direction`` is the local crack propagation direction and ``normal`` is the crack-face normal.
Either vector may be omitted; the missing vector is reconstructed with the ordered crack-front
tangent. If both are omitted, 4C reconstructs the advance direction from the adjacent crack-face
boundary faces. This requires a unique, conforming crack-face geometry at the front.

Contours and local coordinates
------------------------------

For each front node, the local frame is

* :math:`\boldsymbol{x}`: crack propagation direction,
* :math:`\boldsymbol{y}`: crack-face normal,
* :math:`\boldsymbol{z}`: crack-front tangent.

Contour layer 1 contains every solid element incident to a crack-front node. Each later layer adds
elements sharing at least one node with the preceding layer. Therefore, the algorithm works on
unstructured meshes and permits open or truncated domains at free surfaces and symmetry planes.
No geometric radius or closed contour ring is assumed.

Equivalent domain integral
--------------------------

For front node :math:`a` and contour :math:`c`, 4C constructs a scalar virtual-extension field
:math:`q_{a,c}`. Harmonic extension first distributes the crack-front nodal basis into the contour.
The fields are then normalized to form a partition of unity at every contour-interior node, while
remaining zero at the artificial outer contour boundary. Consequently, their sum is one throughout
the interior and drops to zero only across the outer element layer, as required by the standard EDI
construction. Physical free, crack, and symmetry boundaries remain natural boundaries.

The configurational-force vector is evaluated with standard element quadrature as

.. math::

   \boldsymbol{G}_{a,c} =
   \int_{\Omega_c}
   \left((\nabla\boldsymbol{u})^\mathsf{T}\boldsymbol{\sigma}
   - W\boldsymbol{I}\right)\nabla q_{a,c}\,\mathrm{d}V .

It is normalized by the crack-front finite-element support
:math:`A_a=\int_{\Gamma_f}N_a\,\mathrm{d}s`. The reported scalar is

.. math::

   J_{a,c} = \frac{\boldsymbol{G}_{a,c}}{A_a}\cdot\boldsymbol{x}_a .

The predicted maximum-configurational-force direction is the normalized projection of
:math:`\boldsymbol{G}_{a,c}/A_a` onto the plane normal to :math:`\boldsymbol{z}_a`. A vanishing
projected force has no defined propagation direction and is reported as an error.

Output
------

ASCII output creates one appendable YAML file per crack front. Binary output creates one VTK
unstructured-grid time series per crack front, with the front represented by line cells and point
arrays for each selected contour. The structural VTK output additionally contains
``fracture_domain_<crack_id>`` cell data. Its value is zero outside the contour domain and otherwise
the first topological layer containing the element.

``output_contours`` accepts ``all``, ``outer``, or an explicit list of one-based contour numbers.
All contours are evaluated even when only a subset is serialized.

Supported regime and limitations
--------------------------------

The current EDI kernel supports three-dimensional displacement-based solid elements with
``KINEM: linear``. It uses the material's strain-energy density and currently excludes finite
kinematics, EAS and F-bar formulations, shells, prestress, thermal/configurational source terms,
body-force corrections, material interfaces, and non-traction-free crack faces. Unsupported solid
formulations fail explicitly.

``G`` is available only for linear-elastic St. Venant-Kirchhoff materials and is then equal to
``J``. Stress-intensity factors, interaction integrals, and crack growth are not implemented.

Symmetry
--------

4C reports the integral of the modeled domain and never applies an automatic symmetry factor. In
particular, a half model whose crack plane is a symmetry plane is expected to report half the J
value of the corresponding full model. Any desired scaling is the user's responsibility.

M(T) analytical reference
--------------------------

For a center-cracked M(T) specimen with half crack length :math:`a`, half specimen width
:math:`W`, and nominal tensile stress :math:`\sigma`, the plane-strain reference used for
verification is

.. math::

   K_I = \sigma\sqrt{\pi a}
   \sqrt{\sec\left(\frac{\pi a}{2W}\right)}, \qquad
   J = \frac{K_I^2}{E'}, \qquad
   E' = \frac{E}{1-\nu^2}.

For displacement-controlled models, :math:`\sigma` must be calculated from the tensile reaction
resultant divided by the corresponding gross loaded area; it must not be inferred from the imposed
displacement using the uncracked homogeneous strain.

Williams patch verification
----------------------------

The Mode-I disk patch ``crack_tip_disk_williams.4C.yaml`` prescribes the plane-strain Williams
displacement field for :math:`K_I=1000`, :math:`E=210000`, and :math:`\nu=0.3`. The analytical
reference is

.. math::

   J = \frac{K_I^2}{E'} = 4.333333.

The supplied radial mesh uses 32 angular sectors and 10 radial layers. Its collapsed crack-tip
hexahedra are represented as wedges in ``crack_tip_disk_wedge.msh`` because the solid-element
Jacobian must be positive at every element node. Contours 2--8 produce :math:`J=4.370` to
:math:`4.467`; contours 3--8 vary by 0.62 percent and their outer values are 3.1 percent above the
analytical reference. The finite-element displacement has a 2.50 percent relative nodal
:math:`L_2` error, concentrated in the innermost radial layer, and matches the prescribed outer
boundary to roundoff.