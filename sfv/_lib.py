from sfv._sfv import sfvmod

__all__ = [
    "__version_tuple__",
    "predict_spatial_flux",
]

__version_tuple__ = tuple(sfvmod.sfv_v)


def predict_spatial_flux(
    siga0,
    siga1,
    nsf0,
    nsf1,
    keff0,
    adjointMoments,
    forwardMoments,
    lambdaEig,
    forwardFlux0,
    keff1=-1,
    overwrite_flux=False,
):
    """Predict the new spatial flux distribution using SFV

    All arrays are expected to be of type `numpy.float64`.

    Vectors must be ordered such that ``siga0[i]``, ``siga1[i]``,
    ``nsf0[i]``, ``nsf1[i]``, and ``phi[i]``  all correspond to the
    same node. Flux modes ``adjointMoments[i, j]`` and
    ``forwardMoments[i, j]`` represent the ``j``-th moment of the
    adjoint and forward flux, respectively, in node ``i``, similarly
    ordered. ``eig[j]`` is the ``j``-th eigenvalue.

    Vectors must be consistently sized, such that ``N`` nodes and ``M``
    modes of the flux, cross section vectors are ``N``-valued real
    vectors, ``eig`` is a ``M``-valued real vector, and forward and
    adjoint mode matrices are ``NxM`` matrices of reals.

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
        ``k``-eigenvalue at the known point
    adjointMoments : numpy.ndarray
        Adjoint modes ``[i, j]`` of the flux taken from the known
        point. Must be in fortran ordering (column-major)
    forwardMoments : numpy.ndarray
        Forward modes ``[i, j]`` of the flux taken from the known
        point. Must be in fortran ordering (column-major)
    lambdaEig : numpy.ndarray
        Vector of the first ``nModes`` eigenvalues from the
        fission matrix, such that ``lambdaEig[0] == 1/keff0``
    forwardFlux0 : numpy.ndarray
        Fundamental mode of the forward flux, normalized, and
        ordered similarly to ``nsf0`` and other cross sections.
    keff1 : float, optional
        Criticality at next time step. If less than zero,
        then estimate using first order perturbation theory
    overwrite : bool
        If True, then update the phi0 vector in place.
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

    if not overwrite_flux:
        forwardFlux0 = forwardFlux0.copy()

    sfvmod.predict_spatial_flux(
        siga0,
        siga1,
        nsf0,
        nsf1,
        keff0,
        adjointMoments,
        forwardMoments,
        lambdaEig,
        forwardFlux0,
        keff1,
    )
    return forwardFlux0
