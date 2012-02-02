.. highlight:: c++

Distributing local assembly between threads in OP2
==================================================

.. sectionauthor:: Graham Markall <grm08@doc.ic.ac.uk>

Introduction
------------

In this document we consider how we can change the OP2 syntax and semantics in
order to support the automatic distribution of the computations in the local
assembly between threads.

Example problem
...............

Consider the form of the mass matrix:

.. math:: M[i,j] = \int_\Omega \phi_i \phi_j \mathrm{d}X

This defines the value of each entry of the mass matrix. When performing
integration using Gaussian quadrature, this form becomes:

.. math:: M[i,j] = \sum_q w_q \phi_i \phi_j

We might implement the assembly in an OP2 kernel for assembly over triangles
with three quadrature points as: ::

  void mass(double **x, double **mat)
  {
    // Calculate Jacobian
    double J_00 = x[1][0] - x[0][0];
    double J_01 = x[2][0] - x[0][0];
    double J_10 = x[1][1] - x[0][1];
    double J_11 = x[2][1] - x[0][1];

    // Calculate determinant of Jacobian
    double detJ = J_00*J_11 - J_01*J_10;

    // Quadrature weights
    double w[3] = { 0.166667, 0.166667, 0.166667 };

    // Values of basis functions at quadrature points
    double CG1[3][3] = {{ 0.666667, 0.166667, 0.166667 },
                        { 0.166667, 0.666667, 0.166667 },
                        { 0.166667, 0.166667, 0.666667 }};


    // Assembly computation
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
        for (int q=0; q<3; ++q)
          mat[i,j] += CG1[i][q] * CG1[j][q] * detJ * w[q];
  }

This has the problem that it fixes the mapping of elements to threads 1:1.
This is acceptable for low-order polynomials, but higher-order polynomials
will require much more computation per element, and the 1:1 mapping becomes
inefficient. The proposed solution to this problem involves making the kernel
free in the loop indices (``i``, ``j``, and ``q``).

Proposal
--------

Instead of writing the kernel containing the loop, the kernel now contains: ::

  void mass(double **x, double *mat, int i, int j, int q)
  {
    // Calculate Jacobian
    double J_00 = x[1][0] - x[0][0];
    double J_01 = x[2][0] - x[0][0];
    double J_10 = x[1][1] - x[0][1];
    double J_11 = x[2][1] - x[0][1];

    // Calculate determinant of Jacobian
    double detJ = J_00*J_11 - J_01*J_10;

    // Quadrature weights
    double w[3] = { 0.166667, 0.166667, 0.166667 };

    // Values of basis functions at quadrature points
    double CG1[3][3] = {{ 0.666667, 0.166667, 0.166667 },
                        { 0.166667, 0.666667, 0.166667 },
                        { 0.166667, 0.166667, 0.666667 }};


    // Assembly computation
    *mat += CG1[i][q] * CG1[j][q] * detJ * w[q];
  }

The OP2 runtime is now working with the iteration space of the set of
elements, and also iterates over the sub-sets of the basis functions and
quadrature points. In order to inform the OP2 runtime about the iteration
space, we need to modify the ``op_par_loop`` syntax to provide an *iteration
space* rather than just a set. The syntax of the ``op_par_loop`` call in this
case is: ::

  op_par_loop(mass, op_iteration_space(elements, dofmap, dofmap, op_range(0,2)),
               op_arg(x, ... OP_ALL),
               op_arg(mat, op_i[0], op_i[1], OP_INC),
               op_arg(NULL, op_i[2])
             );

Instead of just passing a single set as the second argument, we now pass an
iteration space, implemented by the ``op_iter`` function, which consists of a
set, and then a set of maps that can be used to index into subsets. In the
above example:

* ``elements`` is the set of elements. This first argument is the set that OP2
  will iterate over, as with the usual syntax.

* ``dofmap`` is the mapping from elements to degrees of freedom. From this,
  the OP2 runtime constructs an \emph{indexing object} that ranges over an
  iteration space based on the DOFs that are mapped from a particular element.
  For example, if ``dofmap`` maps element 3 to DOFs 4, 6, and 7, the resulting
  indexing object is a function that maps :math:`0 \rightarrow 4`, :math:`1
  \rightarrow 6`, and :math:`2 \rightarrow 7`. This indexing object is not
  exposed to the user, but is used internally by the OP2 runtime.

* ``op_range(0,2)`` provides an indexing object that provides the identity
  mapping, ranging from 0 to 2.  This construction is necessary to index into
  the array of quadrature points, which don't have any relationship with the
  sets of elements or DOFs, but instead are an arbitrary array inside the
  kernel.

Inside each ``op_arg``, the indexes that are defined in the iteration space
can be referred to by the ``op_i`` array. For example, ``op_i[0]`` tells the
runtime to use the first index based on ``dofmap``, ``op_i[1]`` tells the
runtime to use the second index based on ``dofmap``, and ``op_i[2]`` tells the
runtime to use the index based on ``op_range(0,2)`` - the index into ``op_i``
corresponds to the argument to ``op_iter``, starting at 0 for the second
argument (because the first argument was the set over which the iteration
takes place).

Also, instead of the matrix being passed as a ``double**``, it is passed as a
pointer to a single scalar variable to be updated - a particular term of the
matrix that is being assembled.

Now, as well as the runtime performing the iteration over the elements, it
additionally is in control of the iteration over both the basis function
indexes (given by ``dofmap``, referred to as ``i`` and ``j`` in the kernel)
and also over the quadrature points (referred to as ``q`` in the kernel). As
with the iteration over the set, it is free to iterate over any of these
indices in any order, which also allows OP2 the freedom to map threads to
matrix elements in any way (e.g. one thread per element, one block per
element, or in-between).

Further example
---------------

In this example we also read from another dat - this is an example of the RHS
evaluation. ::

  void rhs(double **x, double *vec)
  {
    // Calculate Jacobian
    double J_00 = x[1][0] - x[0][0];
    double J_01 = x[2][0] - x[0][0];
    double J_10 = x[1][1] - x[0][1];
    double J_11 = x[2][1] - x[0][1];

    // Calculate determinant of Jacobian
    double detJ = J_00*J_11 - J_01*J_10;

    // Quadrature weights
    double w[3] = { 0.166667, 0.166667, 0.166667 };

    // Values of basis functions at quadrature points
    double CG1[3][3] = {{ 0.666667, 0.166667, 0.166667 },
                        { 0.166667, 0.666667, 0.166667 },
                        { 0.166667, 0.166667, 0.666667 }};


    // Assembly computation
    for (int i=0; i<3; ++i)
        for (int q=0; q<3; ++q)
          vec[i] += CG1[i][q] *  * detJ * w[q];
  }

