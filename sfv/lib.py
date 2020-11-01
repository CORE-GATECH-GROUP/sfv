"""
Low-level interface directly to the FORTRAN routines

.. warning::

    These functions provide little to no consistency
    checking, nor error checking. If a safer mode
    of operation is desired, use :func:`sfv.applySFV`


.. autodata:: __version_tuple__
    :annotation: = ``(int, int)`` with the major and minor library version

.. autodata:: __version__
    :annotation: = String representation of the library version

.. autofunc:: predict_spatial_flux

"""


from ._lib import *  # noqa: F403

__version__ = "{}.{}".format(*__version_tuple__)  # noqa: F405
