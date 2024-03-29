.. highlight:: c++

.. _vector_maps:

Vector Maps
===========

.. sectionauthor:: Florian Rathgeber <f.rathgeber10@imperial.ac.uk>

We propose an extension to the OP2 API to allow vector valued maps (i.e. 
maps that associate a single entity in the *source* set to several 
entities in the *target* set as they already exist in OP2) to be passed 
as a single parameter to the kernel and the parallel loop call. This 
affects both the *user  API* at the ``op_par_loop`` invocation and the 
*kernel API* as outlined in this proposal.

Furthermore we propose to eliminate the assumption that data of an 
``op_dat`` associated with a single element of an ``op_set`` is laid out 
contiguously in memory and can hence be addressed by subscripting a 
pointer to a primitive data type. Instead, OP2 should have the freedom 
to choose the data layout and the kernel is provided with a way to 
access the data without knowing the data layout.

.. _vector_maps_motivation:

Motivation
----------

Suppose we want to iterate over elements of a mesh in a parallel loop
and read the coordinates of all the degrees of freedom (DOFs)
associated with each of the elements. We declare a mapping from
*elements* to *DOFs* of the mesh which we call a *vector map* since it
maps a *single source entitiy* (in this case an element) to *a vector
of destination entities* (in this case the DOFs associated with that
element). The *target dimension* (i.e. the size of the vector of target
entities mapped to from a single source entry) is 3 in this example.

Assume OP2 data structures for the example above declared as ::

  op_set p1       = op_decl_set(4, "nodes");
  op_set elements = op_decl_set(2, "elements");

  int p1_ele_i[] = {1,2,3,
                    3,2,4};
  op_map p1_ele = op_decl_map(elements, p1, 3, p1_ele_i, "element-node");

  op_dat coordinates = op_decl_dat(p1,       2, "double", x, "coordinates");
  op_dat com         = op_decl_dat(elements, 2, "double", c, "com");

The call invoking the parallel loop over the set of elements for a
kernel reading ``coordinates`` for each of the 3 DOFs of each element
computing the center of mass ``com`` is given as ::

  op_par_loop ( kernel, "kernel", elements,
                op_arg_dat(coordinates, 0, p1_ele, 2, "double", OP_READ),
                op_arg_dat(coordinates, 1, p1_ele, 2, "double", OP_READ),
                op_arg_dat(coordinates, 2, p1_ele, 2, "double", OP_READ),
                op_arg_dat(com,        -1,  OP_ID, 2, "double", OP_INC)
              );

Arguments related to the mapped data ``coordinates`` need to be
repeated three times for each index into the vector map
``p1_ele`` with the target dimension 3, resulting in 12
parameters for the mapped data. Analogously, for hexahedrons instead of
triangles the vector map is of target dimension 6 such that the number
of parameters for the mapped data increases to 24. It is obvious that
the number of parameters becomes unmanageable for vector maps of large
target dimension that would arise e.g. in high-order finite element
methods. Furthermore, since it is usually the case that all data items
targeted by a vector map are accessed in a parallel loop and all in the
same way, it is counter-intuitive to have to give the mapping multiple
times.

The same argument applies to the kernel, which gets passed a pointer
for each component in the target dimension of the vector map ::

  void kernel(const float *x0, const float *x1, const float *x2, float *c) {
    c[0] = (x0[0] + x1[0] + x2[0])/3.0f;
    c[1] = (x0[1] + x1[1] + x2[1])/3.0f;
  }

Specification
-------------

In the previous example we associated a 2D coordinate with every degree
of freedom of the mesh, i.e. an ``op_dat`` with a two component vector
for every element of the associated ``op_set``. Similarly we could
imagine higher rank tensor-valued data such as a diffusivity tensor
associated with each element of an ``op_set``.

Data to the kernel given in the previous example is passed as a pointer
to a primitive type, assumed to point to a memory location where it is
stored contiguously. However, optimal storage layouts of this data
depend on the target architecture and there is not a single best
solution that fits all cases. Consequently, OP2 needs to have the
freedom of choosing a data layout tailored to the architecture and
abstract this choice away from the user by providing uniform access
methods.

Requirements
~~~~~~~~~~~~

- There is a facility to select all components in the target dimension
  of a vector map when passing indirectly accessed data to a parallel
  loop. The number of parameters passed to a parallel loop call and to
  a kernel relating to a vector map are independent of the target
  dimension. Furthermore, there is no difference in the number of
  parameters for specifying directly accessed data and data accessed
  through a vector map.
- Data associated with a single element of a set may be
  multi-dimensional. No assumptions about the data layout can be made,
  OP2 has the implementation choice.
- Data access functions uniformly abstract the choice of data layout.
  A scalar item of data can be located by a tuple of indices. There is
  an index selecting the component of the vector map (if the map is
  multi-dimensional). There is a number of indices equal to the rank of
  data associated with each member of the set (none if it is scalar).
- The target dimension of a vector map and the shape of data associated
  with a single element of a set are known compile-time constants.
- Compatibility to the legacy API is preserved.
- The kernel :term:`ABI` is interoperable between C++ and Fortran. This
  allows any kernel written in either language to be called from OP2
  without the need for a distinction.

The User API
~~~~~~~~~~~~

Recall the example in :ref:`vector_maps_motivation`, where a vector map 
of target dimension 3 was declared, mapping each element to its three 
associated degrees of freedom. This required passing the data, the 
mapping and the access descriptor three times for each index 0, 1, and 2 
into the vector map.

We propose introducing an index specifier ``OP_ALL`` to signal that all 
components of a vector map are accessed in a parallel loop. The ``OP_ALL`` 
specifier serves the role of a vector index similar in semantic to the 
``:`` in Fortran. This reduces the number of parameters to four for each 
``op_dat`` accessed in the parallel loop, irrespective of whether it is 
accessed directly or via a mapping and what the target dimension of that 
mapping is.

The parallel loop call from the previous example can thereby be written 
as follows ::

  op_par_loop ( kernel, "kernel", elements,
                op_arg_dat(coordinates, OP_ALL, p1_ele, 2, "double", OP_READ),
                op_arg_dat(com,         -1,     OP_ID,  2, "double", OP_INC)
              );

This extension does not eliminate the option to specify only a certain
index into a vector map in the parallel loop call and is hence
compatible with the legacy API. Furthermore it does not require any
changes to the way mappings are currently implemented in OP2.

The Kernel API
~~~~~~~~~~~~~~

OP2 kernels are written against a public kernel API that is compiled 
unaltered for the serial reference implementation only. For parallel 
execution on a target platform, the kernel is processed with a 
source-to-source translator and the kernel body embedded in 
platform-specific kernel that takes care of data marshaling etc. as 
necessary. Performance of the public kernel API implementation is hence 
not important and we can safely require data passed to the kernel to be 
in a suitable dense storage layout populated by a copy-in and written 
back by copy-out if necessary.

The proposal requires the kernel API to be changed in two respects. 
Indirectly accessed data is passed to the kernel as one parameter 
instead of one per target dimension of the vector map. This parameter 
must furthermore be of a type that is subscriptable with a number of 
indices according to the rank of the data as described in the 
requirements.

Maintaining ABI compatibility between C++ and Fortran
.....................................................

We can achieve ABI compatibility between C++ and Fortran by passing 
vector-mapped data to OP2 kernels as:

 - a multi-dimensional array (dope vector) in Fortran, and
 - a fixed-size multi-dimensional array in C.

This amounts to effectively passing a pointer to a primitive type (e.g. 
a ``double*``) to the kernel in both cases and hence is ABI interoperable 
between C++ and Fortran and also compatible to the existing kernel ABI 
for the cases of directly accessed data.

Kernel API implementation
.........................

In the kernel body, data can be addressed with a number of indices equal
to the rank of the mapped data plus an index selecting the component of
the vector map (if applicable). The actual data pointed to is in both
cases assumed to be stored densely in row-major array format in C and
column-major array format in Fortran. It is assumed to be populated by
copy-in before kernel invocation and written back by copy-out afterward
if necessary. 

Multi-dimensional arrays are a built-in feature of the Fortran language. 
The example kernel using this data abstraction is given below:

.. code-block:: fortran

  subroutine kernel ( x , c ) BIND(C)
    real ( 8 ) , dimension (3,2) :: x
    real ( 8 ) , dimension (2) :: c
    c(0) = (x(0,0) + x(1,0) + x(2,0))/3.0
    c(1) = (x(0,1) + x(1,1) + x(2,1))/3.0
  end subroutine kernel

C supports multi-dimensional arrays only if all the bounds (the shape)
is known at compile time. This is the case for OP2 kernel parameters and
hence not a restriction. The kernel from the example using this data
abstraction can be written as: ::

  void kernel( const double x[3][2], double c[2]) {
    c[0] = (x[0][0] + x[1][0] + x[2][0])/3.0f;
    c[1] = (x[0][1] + x[1][1] + x[2][1])/3.0f;
  }

Complete example using the specified API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following listing shows a complete example using the specificied C++
API.  It contains a kernel and an OP2 main program with a parallel loop
call invoking the kernel.

.. literalinclude:: vectormaps_example.cpp
