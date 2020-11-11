sfv
===
Small library for predicting variations in the scalar neutron flux using perturbation theory

.. image:: https://zenodo.org/badge/DOI/10.1080/00295639.2019.1661171.svg
   :target: https://doi.org/10.1080/00295639.2019.1661171
   :alt: Nuclear Science and Engineering 10.1080/00295639.2019.1661171

This repository holds the code necessary to build the ``sfv`` Python package.
The functionality is quite small, providing three simple functions:

* ``sfv.applySFV``: Apply the prediction to obtain normalized scalar flux
* ``sfv.getAdjFwdEig``: Extract ``k`` modes of the foward and adjoint fission
  source and associated $k$ eigenvalues from the fission matrix
* ``sfv.lib.predict_spatial_flux``: Direct access to the SFV compiled library with
  a similar API as ``sfv.applySFV``, but with little to no safety checks.

Installation
------------

Installing ``sfv`` requires ``numpy`` and ``f2py`` in order to compile the FORTRAN
library. The preferred way to install is using a virtual environment, as this produces
a system that is easier to upgrade in the future. Using the ``venv`` module, the steps
to create a virtual environment are

.. code-block:: shell

   $ python -m venv /path/to/venv
   $ source /path/to/venv/bin/activate
   (venv) $ which python
   /path/to/venv/bin/python

Note that python >= 3.5 is supported for this project.

Dependencies can be installed with

.. code-block:: shell

   (venv) $ pip install -r requirements.txt

and the entire package can be built and installed with

.. code-block:: shell

   (venv) $ pip install .

Alternatively, the following two commands will build the compiled
library into the ``build`` directory and then install the package

.. code-block:: shell

   (venv) $ python setup.py build
   (venv) $ python setup.py install

During the build stage, ``numpy`` and the ``f2py`` backend will attempt to find existing
FORTRAN compilers, as well as the Intel Math Kernel Library (MKL) if it exists. The
implementation is not so demanding that MKL is required, but it could provide some benefits.

Upgrading / uninstalling
~~~~~~~~~~~~~~~~~~~~~~~~

Unforuntately, the linking of the compiled library to the python library requires
``distutils`` rather than ``setuptools`` meaning that the entire package must be deleted
by hand when re-installing or upgrading. For further discussion, see
`this issue in the PYPA GitHub <https://github.com/pypa/pip/issues/5247#issuecomment-381550610>`_

Testing
-------

Test coverage is not currently great here in this repository. However, the code has been
tested on practical problems with good results. Eventually, these data should be added
to this repository to ensure the robustness of the package.

Testing is handled with ``pytest``, which is the only dependency listed
in ``requirements-test.txt``. It can then be installed with

.. code-block:: shell

   (venv) $ pip install pytest

or

.. code-block:: shell

   (venv) $ pip install -r requirements-test.txt

Tests should be run from the ``tests`` directory

.. code-block:: shell

   (venv) $ cd tests
   (venv) $ pytest

since the python interpreter will look for the compiled library in the ``sfv``
directory, not in the virtual environment.

FORTAN compatibility
--------------------

The implementation of the SFV prediction is written in a FORTRAN-90 file contained in
``lib/sfv.f90``. For transparency, the library was written to interface with python
first and (maybe) FORTRAN second. If you have issues with linking this library in a
FORTRAN program, or have advice on this, we welcome that insight.

Citing
------

If you use this package, please cite the above *Nuclear Science and Engineering*

.. code-block:: bibtex

   @Article{<your bibtex key here>
      author  = {Johnson, Andrew E. and Kotlyar, Dan},
      title   = {A Transport-Free Method for Predicting the Post-Depletion Spatial Neutron Flux Distribution},
      doi     = {10.1080/00295639.2019.1661171},
      number  = {2},
      pages   = {120-137},
      url     = {https://doi.org/10.1080/00295639.2019.1661171},
      volume  = {194},
      journal = {Nuclear Science and Engineering},
      year    = {2020},
   }

If you'd like, let us know about your publication and we'll add it
to a (to be created) publication list.