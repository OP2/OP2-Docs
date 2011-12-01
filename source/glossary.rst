.. _glossary:

Glossary
========

.. glossary::

  ABI (Application Binary Interface)
    The low-level interface (i.e. the type of binary data passed on a function call) of a library or software program.

  API (Application Programming Interface)
    The public interface (i.e. set of function signatures) exposed by a library or software program the API user needs to follow.

  extent
    The size of a given :term:`rank` of a multi-dimensional array.

  OP2 global API
    The API used to write an OP2 program, providing declaration of OP2 data types (set, data, map) and parallel loop calls.

  OP2 kernel API
    The API defining the signature of an OP2 kernel. Kernels using this API will be compiled unmodified (without undergoing source-to-source translation) only for the serial reference case.

  rank
    The number of subscripts to index into a multi-dimensional array.

  shape
    The tuple of :term:`extent` s corresponding to each :term:`rank` of a multi-dimensional array giving its overall size.

  target dimension
    The size of the vector of target entities mapped to from a single source entry of a :term:`vector map`.

  vector map
    A mapping from a *single source entitiy* to *a vector of destination entities*. We call the size of the vector of target entities mapped to from a single source entry the :term:`target dimension`.
