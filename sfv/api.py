from numbers import Integral, Real

import numpy
from numpy.linalg import LinAlgError
from scipy.sparse.linalg import eigs as speig
from scipy.sparse import issparse, csr_matrix, spmatrix

from .lib import (
    __version_tuple__,
    predict_spatial_flux,
)

__all__ = [
    "getAdjFwdEig",
    "applySFV",
]

_expected = (0, 3)

# compare linked fortran module to expected
if __version_tuple__[0] != _expected[0] or __version_tuple__[1] < _expected[1]:
    raise ImportError(
        "Linked sfv fortran library {} is incompatible with expected "
        "version {}".format(__version_tuple__, _expected)
    )
del _expected


def applySFV(
    siga0,
    siga1,
    nsf0,
    nsf1,
    keff0,
    adjMom,
    fwdMom,
    eig,
    phi0,
    keff1=None,
    overwrite=False,
):
    """
    Apply the spatial flux variation method to predict normalized flux

    ``siga0`` and ``nsf0`` correspond to the known time point,
    typicaly beginning-of-step, where ``keff0``, ``adjMom``,
    ``fwdMom``, and ``eig`` were calculated. ``siga1``, and ``nsf1``
    correspond to the later point in time, where the flux prediction
    will occur.

    Vectors must be ordered such that ``siga0[i]``, ``siga1[i]``,
    ``nsf0[i]``, ``nsf1[i]``, and ``phi[i]``  all correspond to the
    same node. Flux modes ``adjMom[i, j]`` and ``fwdMom[i, j]``
    represent the ``j``-th moment of the adjoint and forward flux,
    respectively, in node ``i``, similarly ordered. ``eig[j]``
    is the ``j``-th eigenvalue. Vectors must be consistently sized,
    such that ``N`` nodes and ``M`` modes of the flux, cross
    section vectors are ``N``-valued real vectors, ``eig`` is a
    ``M``-valued real vector, and forward and adjoint mode matrices
    are ``NxM`` matrices of reals.

    Forward and adjoint matrices, as well as eigenvalues, will be
    coerced to real valued quantities. If reals are given, the
    imaginary part may be dropped.

    Parameters
    ----------
    siga0 : numpy.ndarray
        Macroscopic absorption cross section at known point
    siga1 : numpy.ndarray
        Macroscopic absorption cross section at prediction point
    nsf0 : numpy.ndarray
        Macroscopic nu-sigma fission cross section at known point
    nsf1 : numpy.ndarray
        Macroscopic nu-sigma fission cross section at prediction point
    keff0 : float
        Criticality at known point
    adjMom : numpy.ndarray
        Adjoint modes ``[i, j]`` of the flux taken from the known
        point
    fwdMom : numpy.ndarray
        Forward modes ``[i, j]`` of the flux taken from the known
        point
    eig : numpy.ndarray
        Vector of the first ``nModes`` eigenvalues from the
        fission matrix
    phi0 : numpy.ndarray
        Fundamental mode of the forward flux, normalized, and
        ordered similarly to ``nsf0`` and other cross sections.
    keff1 : None or float
        Criticality at next time step. If None or less than zero,
        then estimate using first order perturbation theory
    overwrite : bool
        If trueish, then update the phi0 vector in place.
        Otherwise, create a copy. The predicted flux
        will still be returned, regardless of this value

    Returns
    -------
    numpy.ndarray
        Normalized flux prediction due to the changes between
        known and predicted states. Ordered identically to the
        macroscopic cross sections, e.g. ``predFlux[j]`` is the
        predicted, normalized flux in the same node as
        ``siga1[j]``.

    """

    # Sanitize and check for consistency

    siga0, siga1, nsf0, nsf1 = _prepSfvXs(siga0, siga1, nsf0, nsf1)

    eig = numpy.asarray(eig, dtype=numpy.float64)
    if len(eig.shape) != 1:
        eig = eig.squeeze()

    adjMom = numpy.asfortranarray(adjMom, dtype=numpy.float64)
    if adjMom.shape != (siga0.size, eig.size):
        raise ValueError(
            "Adjoint flux moment matrix not consistent with XS nor "
            "eigenvalues. Flux shape {}, cross section {}, number of "
            "eigenvalues {}".format(adjMom.shape, siga0.size, eig.size)
        )

    fwdMom = numpy.asfortranarray(fwdMom, dtype=adjMom.dtype)

    if fwdMom.shape != adjMom.shape:
        raise ValueError(
            "Flux moments are inconsistent. Adjoint : {}. "
            "Forward: {}".format(adjMom.shape, fwdMom.shape)
        )

    if keff1 is None:
        keff1 = -1.0
    elif not isinstance(keff1, Real):
        raise TypeError("Keff - {}".format(type(keff1)))

    phi0 = numpy.asarray(phi0, dtype=numpy.float64)
    if phi0.shape != siga0.shape:
        raise ValueError(
            "Forward flux inconsistent with cross sections. Expected {}, "
            "got {}".format(siga0.shape, phi0.shape)
        )

    return predict_spatial_flux(
        siga0,
        siga1,
        nsf0,
        nsf1,
        keff0,
        adjMom,
        fwdMom,
        eig,
        phi0,
        keff1=keff1,
        overwrite_flux=overwrite,
    )


def _prepSfvXs(*arrays):
    out = tuple(numpy.asarray(x, dtype=numpy.float64) for x in arrays)
    shapes = {x.shape for x in out}
    if len(shapes) > 1:
        raise ValueError(
            "SFV XS vectors are not of consistent shape. " "Found {}".format(shapes)
        )

    shape = shapes.pop()
    if len(shape) != 1:
        raise ValueError("SFV XS vectors are not vectors, instead are {}".format(shape))
    elif not shape[0]:
        raise ValueError("SFV XS vectors are empty")

    return out


def getAdjFwdEig(A, numModes=None):
    """
    Compute adjoint and foward fission source moments, and eigenvalues

    Parameters
    ----------
    A : scipy.sparse.spmatrix
        Double precision square fission matrix. Shape should be ``NxN``
    numModes : int, optional
        Number of modes to be extracted. If given, must be less than
        ``N-1``

    Returns
    -------
    adj : numpy.ndarray
        Adjoint moments from of the fission source
    fwd : numpy.ndarray
        Forward moments from of the fission source
    eig : numpy.ndarray
        $k$-eigenvalues of the fission matrix

    Raises
    ------
    numpy.linalg.LinAlgError
        QR algorithm for the eigensolver failed and no eigenvectors
        were computed

    """
    if not issparse(A):
        A = csr_matrix(A)
    if not A.ndim == 2:
        raise ValueError(f"A must be 2D square matrix, not {A.shape}")
    (nr, nc) = A.shape
    if nr != nc:
        raise ValueError(f"A must be 2D square matrix, not {A.shape}")

    if numModes is None:
        numModes = nr - 2
    elif not isinstance(numModes, Integral):
        raise TypeError(f"Number of modes k must be positive integer, not {numModes}")
    elif numModes <= 0:
        raise ValueError(f"Number of modes k must be positive integer, not {numModes}")

    try:
        fwdMoments, fwdKEigs = _eigWrapper(A, numModes)
    except Exception as ee:
        raise LinAlgError("Failed to obtain forward moments") from ee

    try:
        adjMoments, _adjKEigs = _eigWrapper(A.T, numModes)
    except Exception as ee:
        raise LinAlgError("Failed to obtain adjoint moments") from ee

    return adjMoments, fwdMoments, fwdKEigs


def _eigWrapper(A: spmatrix, k: int):
    # Take magnitude and sort according to largest k eigenvalues
    keigs, vectors = speig(A, k, which="LM", return_eigenvectors=True)
    keigs = keigs.real
    ix = keigs.argsort()[::-1]
    return vectors[:, ix].real, keigs[ix]
