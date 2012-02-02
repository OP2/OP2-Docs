.. highlight:: c++

Finite Element Computations
===========================

.. sectionauthor:: David A. Ham <david.ham@imperial.ac.uk>

The finite element method primarily consists of integrals over cells and, in
some cases, facets [#]_. In this regard, it looks like a prime contender for
OP2. The complication arises because the degrees of freedom in a finite
element mesh do not all coincide with the same geometric entities. There are
two important cases of this. The first is higher order elements, the second is
discontinuous elements. Figure :ref:`fig:elements` illustrates three discrete
function spaces frequently employed in finite element. These three embody most
of the cases we need to consider in the first instance. In the first instance,
we shall consider explicit methods. Implicit methods will follow as a direct
result of applying sparse matrix functionality to finite element kernels.

.. _fig:elements:

.. figure:: elements.*
   :align: center

   Finite element degrees of freedom

   Finite element degree of freedom locations for a small mesh
   region. The discretisations correspond to piecewise linear, continuous
   (:math:`P1`, left), piecewise quadratic, continuous (:math:`P2`, centre) and
   piecewise linear, discontinuous (:math:`P1_\mathrm{DG}`, right) function spaces
   respectively. The two elements in each figure are labelled :math:`e_1` and
   :math:`e_2` respectively.

Element maps
------------

Let's declare OP2 sets which correspond to the meshes in figure
:ref:`fig:elements`: ::

  op_set p1, p2, p1dg;

  op_set p1(4,NULL,"P1");
  op_set p2(9,NULL,"P2");
  op_set p1dg(6,NULL,"P1DG");

Of course, we're also going to need a set of elements: ::

   op_set elements(2, NULL);

The primary access pattern in the finite element method is via elements, so
we will need maps from the finite element sets to the elements: ::

  int p1_ele_i[] = {1,2,3,
                    3,2,4};
  op_map p1_ele(elements, p1, 3, p1_ele_i,"p1_ele_map");

  int p2_ele_i[] = {1,3,2,6,5,4,
                    4,5,6,7,8,9};
  op_map p2_ele(elements, p2, 6, p2_ele_i,"p2_ele_map");

  int p1dg_ele_i[] = {1,2,3,
                      4,5,6};
  op_map p1dg_ele(elements, p1dg, 3, p1dg_ele_i,"p1dg_ele_map");

The maps need to be internally consistent according to some algorithm,
so that the element has a consistent orientation in all maps. The
orientation here is the Fluidity `One True Element Numbering
<http://amcg.ese.ic.ac.uk/index.php?title=Local:One_True_Element_Numbering>`_
but others are possible. 

There is therefore an assumption that OP2 will preserve the local ordering
of the entries in a map through any renumberings it conducts. 

Element integrals
-----------------

The basic operation in finite element is the computation of integrals over
elements. Let us take as an example, the simplest integral known to finite
element:

.. math:: \int_\Omega v f \mathrm{d} x, \forall v \in V.
   :label: form

In this integral, :math:`\Omega` is the whole domain (ie the whole mesh),
:math:`f` is some function defined over some discrete function space
:math:`F`. :math:`V` is a basis for a discrete function space (which may or
may not be the same space as :math:`F`). Finite element function spaces are
finite-dimensional and each has a known basis. For current purposes, we can
identify a basis with a set of degree of freedom values. This means that
equation :eq:`form` consists of one integral for each degree of freedom in the 

The usual way to construct the integrals is to decompose :math:`\Omega` into a
set of elements. We can then construct the ``op_dat`` corresponding to ``V``
by looping over each element and computing the contribution of the values of
:math:`f` in that element to the integrals of those ``v`` which touch that
element. This would appear as the following kernel invocation: ::

  op_par_loop(integral_kernel, elements,
              v_dat, 0, v_ele_map, OP_INC, &
              f_dat, 0, f_ele_map, OP_READ, &
              x_dat, 0, x_ele_map, OP_READ);

Were the ``dat`` and ``ele_map`` arguments have been associated with the
data and element maps for :math:`V`, :math:`F` and :math:`X` respectively.

Assembling a matrix from an element integral
--------------------------------------------

The simplest matrix case we can consider is the mass matrix. This matrix is
defined such that:

.. math::  \mathrm{M}f = \int_\Omega v f \mathrm{d} x, \forall v \in V.

The mass matrix assembly would appear in OP2 as: ::

  op_sparsity mass_sparsity(v_ele_map, f_ele_map);
  op_mat<double> mass_mat(op_sparsity);

  op_par_loop(mass_kernel, elements,
              mass_mat, 0, v_ele_map, 0, f_ele_map, OP_INC, &
              x_dat, 0, x_ele_map, OP_READ);

Facet maps
----------

Facet integrals in finite element are similar to face integrals, however
there is some added complexity. There are two maps which are applicable, the
first is for those degrees of freedom whose basis functions are non-zero
on the facet. For current purposes, we can identify those basis functions
with the degrees of freedom lying on the facet. Figure :ref:`fig:facets`
illustrates this arrangement.

.. _fig:facets:

.. figure:: facet.*
   :align: center

   Facet degrees of freedom

   For each mesh, the degrees of freedom in blue contribute to
   simple facet integrals over this facet.

If we consider the discontinuous function space on the right, we can see
that there will be 4 degrees of freedom associated with each interior facet
and two with each exterior facet. This would result in a map with variable
dimension, which is _verboten_ in OP2. However, we can adopt the
FEniCS solution to this problem and separate the surface face integrals from
the interior ones. This also reflects the fact that finite element
techniques usually have to do something different on the boundary. This
means that for our trivial mesh, we'd define: ::

  op_set facets(1,NULL);
  op_set boundary_facets(4,NULL);

  int p1_facet_i[] = {2,3};
  op_map p1_facet_int(facets, p1, 2, p1_facet_i, "p1_facet_map");
  int p1_boundary_i[] = {1,3,
                         3,4,
                         4,2,
                         2,1};
  op_map p1_facet_bnd(boundary_facets, p1, 2, p1_boundary_i, "p1_boundary_map");

  int p2_facet_i[] = {6,5,4};
  op_map p2_facet_int(facets, p2, 3, p2_facet_i, "p2_facet_map");
  int p2_boundary_i[] = {1,2,4,
                         4,7,9,
                         9,8,6,
                         6,3,1};
  op_map p1_facet_bnd(boundary_facets, p1, 3, p2_boundary_i, "p2_boundary_map");
  
  int p1dg_facet_i[] = {2,3,5,4};
  op_map p1dg_facet_int(facets, p2, 4, p1dg_facet_i, "p1dg_facet_map");
  int p1_boundary_i[] = {1,3,
                         4,6,
                         6,5,
                         2,1};
  op_map p1_facet_bnd(boundary_facets, p1, 2, p1_boundary_i, "p1_boundary_map");

.. To do: explain gradient maps.

.. rubric:: Footnotes

.. [#] a facet is a geometric entity of co-dimension 1 so the facets of 3D cells are faces and the facets of 2D cells are edges

