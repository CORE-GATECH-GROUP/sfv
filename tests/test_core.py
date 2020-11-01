import pathlib
from collections import namedtuple

import numpy
import pytest
import sfv


SfvVectorBundle = namedtuple("SfvVectorBundle", "siga0 siga1 nsf0 nsf1 phi0 phi1")
SfvFluxMomentBundle = namedtuple("SfvFluxMomentBundle", "adjoint forward eigenvalues")


class SfvHarness:
    vectorfile = "vectors.txt"
    eigfile = "eigenvalues.txt"
    forwardfile = "forward_moments.txt"
    adjointfile = "adjoint_moments.txt"
    referenceFile = "prediction_reference.txt"
    failFluxFile = "prediction_failure.txt"

    def __init__(self, datadir, keff):
        self.datadir = pathlib.Path(datadir)
        self.keff = keff
        self.xsdata = None
        self.fluxdata = None

    def saveXSData(self, siga0, siga1, nsf0, nsf1, phi0, phi1):
        stack = numpy.stack((siga0, siga1, nsf0, nsf1, phi0, phi1))
        header = "{} {}\n Rows: ".format(*stack.shape)
        header += " ".join(["SIGA0", "SIGA1", "NSF0", "NSF1", "PHI0", "PHI1"])
        numpy.savetxt(self.datadir / self.vectorfile, stack, header=header)
        self.xsdata = SfvVectorBundle(siga0, siga1, nsf0, nsf1, phi0, phi1)

    def loadXSData(self):
        if self.xsdata is None:
            data = numpy.loadtxt(self.datadir / self.vectorfile)
            self.xsdata = SfvVectorBundle(*data)
        return self.xsdata

    def saveFluxMoments(self, fwd, adj, eigv):
        mtxHeader = "{} {}".format(*fwd.shape)
        numpy.savetxt(self.datadir / self.forwardfile, fwd, header=mtxHeader)
        numpy.savetxt(self.datadir / self.adjointfile, adj, header=mtxHeader)
        numpy.savetxt(self.datadir / self.eigfile, eigv.T, header=str(eigv.size))
        self.fluxdata = SfvFluxMomentBundle(adj, fwd, eigv)

    def loadFluxMoments(self):
        if self.fluxdata is None:
            eigenvalues = numpy.loadtxt(self.datadir / self.eigfile)
            forward = numpy.loadtxt(self.datadir / self.forwardfile)
            adjoint = numpy.loadtxt(self.datadir / self.adjointfile)
            assert forward.shape == adjoint.shape
            self.fluxdata = SfvFluxMomentBundle(adjoint, forward, eigenvalues)
        return self.fluxdata

    def runSfv(self, inplace=False, library=False):
        xsdata = self.loadXSData()
        fluxdata = self.loadFluxMoments()

        if library:
            sfvfunc = sfv.lib.predict_spatial_flux
            kws = {"overwrite_flux": inplace}
        else:
            sfvfunc = sfv.applySFV
            kws = {"overwrite": inplace}

        phi0 = xsdata.phi0.copy() if inplace else xsdata.phi0

        prediction = sfvfunc(
            xsdata.siga0,
            xsdata.siga1,
            xsdata.nsf0,
            xsdata.nsf1,
            self.keff,
            fluxdata.adjoint,
            fluxdata.forward,
            1 / fluxdata.eigenvalues,
            phi0,
            **kws,
        )

        if inplace:
            assert prediction is phi0

        return self.comparePrediction(xsdata.phi1, prediction)

    def comparePrediction(self, trueFlux, predFlux):
        refPred = numpy.loadtxt(
            self.datadir / self.referenceFile, unpack=True, skiprows=2
        )[0]
        try:
            assert predFlux == pytest.approx(refPred)
        except AssertionError:
            data = numpy.empty((predFlux.size, 3), dtype=predFlux.dtype, order="F")
            data[:, 0] = predFlux
            numpy.subtract(predFlux, trueFlux, out=data[:, 1])
            rel = data[:, 1].copy()
            nz = trueFlux.nonzero()
            rel[nz] /= trueFlux[nz]
            data[:, 2] = rel
            header = (
                "{} {}\nColumns: Flux, absolute difference, relative difference".format(
                    *data.shape
                )
            )
            numpy.savetxt(self.datadir / self.failFluxFile, data, header=header, fmt="% .18e")
            raise


@pytest.fixture(scope="module")
def sfv10Node():
    harness = SfvHarness(pathlib.Path(__file__).parent / "pin_10node", keff=1.21366)
    harness.loadXSData()
    harness.loadFluxMoments()
    return harness


@pytest.mark.parametrize("inplace", [False, True])
@pytest.mark.parametrize("library", [False, True])
def test_10nodePin(sfv10Node, inplace, library):
    sfv10Node.runSfv(inplace, library)
