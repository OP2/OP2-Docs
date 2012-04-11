This is a repository of specifications for the `OP2 project`_. Eventually it
is supposed to contain all of OP2's documentation. Meanwhile you can find most
of the documentation on the `project page`_ and the `runtime library
repository`_.

To build it you need Sphinx_ (tested with version 1.0.7 and above) and GNU make.

Run

* ``make html`` to build the html docs
* ``make latexpdf`` to build the PDF manual
* ``make publish`` to build both and publish in ``~/public_html/op2``
* ``make publish PUBLISHDIR=<dir>`` to build both and publish in ``<dir>``

.. _project page:
.. _OP2 project: http://people.maths.ox.ac.uk/gilesm/op2
.. _runtime library repository: https://github.com/OP2/OP2-Common
.. _Sphinx: http://sphinx.pocoo.org
