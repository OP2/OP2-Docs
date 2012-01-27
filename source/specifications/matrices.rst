.. highlight:: c++

Matrices
========

In an OP2 context, a (sparse) matrix is a linear operator from one set to
another. In other words, it is a linear function which takes a dataset on
one set ``A`` and returns the value of a dataset on another set ``B``. Of
course, in particular, ``A`` may be the same set as ``B``.

This makes the operation of at least some matrices equivalent to the
operation of a particular OP2 kernel. As an illustrative example, let us
take the example of a finite element operator acting on a set defined on the
vertices of a mesh. Using the local matrix approach described in
[Markall2010]_, we can describe the action of a finite element operator
by the matrix operation:

.. math:: A^T\mathrm{M}^eA.

In this equation, :math:`\mathrm{M}^e` is a block diagonal matrix in which
each diagonal block represents the contribution to the operator from one
element.  The matrix :math:`A` encodes the mapping which selects the nodes for
each element from the set. :math:`A^T` has the effect of mapping the
contributions from each element back to the set.

We can now imagine an OP2 kernel which implements the same finite element
operator. It will be invoked with a call like: ::

  op_par_loop(kernel, elements,
                to_data, map1, index1, OP_INC,
                from_data, map1, index1, OP_READ,
                ..., ..., ..., ..., ..., ...,
                argN, mapN, indexN, accessN)

Here, the kernel plays the role of :math:`\mathrm{M}^e` with the third and
subsequent arguments defining the entries of :math:`\mathrm{M}^e`. Note that
the kernel will be executed once for each element. 

``map1`` plays the role of :math:`A`, it maps from elements to vertices
and is therefore used to select the vertices which participate in the
current element. 

It is therefore possible to see a sparse matrix as a partial evaluation of
this kernel: the kernel is evaluated for the third and subsequent arguments
yielding an operator with two arguments. Since we restrict the kernel to
those which are linear in their first two arguments, we ensure that this
partial evaluation can be described by a set of :math:`\mathrm{dim1}\times
\mathrm{dim1}\times\mathrm{elements}` floats. Indeed, the local matrix
approach consists of storing exactly this set of numbers. To place these
approaches in context, evaluating the full kernel corresponds to a
matrix-free approach: nothing is cached, the full local assembly is
conducted every time the matrix is applied to a set. The local matrix
approach consists of caching the kernel results and applying the maps each
time the matrix is applied to a set. The CSR storage pattern
and other traditional sparse matrix storage schemes consist of taking the
partially evaluated kernel responses and mapping them back onto a direct
linear mapping between the sets. 

It is important to note that nothing in this example is specific to finite
element. All implicit PDE techniques generate sparse matrices of the form
discussed here. Finite element differs from them only in the contents of the
kernels.

Including sparse matrices in OP2
--------------------------------

We can now see that the role of the kernel in constructing our matrix is to
take the third and subsequent arguments and construct the
:math:`\mathrm{dim1}\times \mathrm{dim1}` local tensor. The mapping ``map1``
tells OP2 what the contribution of the local tensor is and it is an
implementation choice what OP2 chooses to do with this information. We might
therefore include a new function which, rather than returning a local vector
through the map, returns a local tensor: ::

  op_par_loop(kernel, elements,
                matrix, map1, index1, map1, index2, OP_INC,
                ..., ..., ..., ..., ..., ...,
                argN, mapN, indexN, typeN, dimN, accessN)

Note that now we have no from and to data, the resulting local tensor is to be
inserted into the opaque matrix object ``matrix``. Suppose that the kernel
returns ``lmatrix``, then for element ``e``, one possible implementation of
the global assembly would be: ::

  call addto(matrix, map1(e), map1(e), lmatrix)

This form of call facilitates the use of third party sparse matrix libraries
as it is of the form which will be expected by them.

Note that matrices operating between two non-identical sets can also be
accommodated in the above code simply by using two different maps.

It is anticipated that the ``matrix`` object will be a wrapper around a third
party library such as Petsc, but nothing prevents the future development of an
OP2 native sparse matrix format if this is thought to be useful. 

Matrix sparsity
---------------

It is a frequent requirement of sparse matrix libraries that the sparsity
pattern be supplied in advance of matrix population. This can be achieved by
a call to a sparsity generation routine: ::

  op_sparsity mat_sparsity(elements, map1, map1);

This will enable the actual matrix object to be instansiated: ::

  op_sparse_matrix<double> matrix(mat_sparsity);

It is likely that the underlying library will require more information than
this, but not yet clear to me what information in particular.

Additional complications
------------------------

Matrices constructed by more than one local assembly operation, for example DG
matrices which have element integrals and facet integerals are a bit more
complicated. The ``op_par_loop`` calls are the same: a matrix can be updated
by more than one loop. For example one loop might operate over faces and
another loop myight operate over elements. This is no different from what may
occur in the construction of a dataset. However, it might be necessary to
extend the ``op_sparsity`` class. For example, this might be expanded to take
more maps. In our DG example, the call might become: ::

  op_sparsity mat_sparsity( &
       elements, vertex_element_map1, vertex_element_map1, &
       edges,    vertex_edge_map1,    vertex_edge_map1)

.. [Markall2010] Markall, Graham R., David A. Ham, and Paul H.J. Kelly.
   "Towards generating optimised finite element solvers for GPUs from high-level specifications."
   Procedia Computer Science 1, no. 1 (May 2010): 1809-1817. 

