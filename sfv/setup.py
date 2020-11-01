def getlibraries():
    """Use system information to find lapack, blas libraries"""
    from distutils.errors import LinkError
    from numpy.distutils import system_info

    hasmkl = system_info.get_info("mkl")
    if hasmkl:
        return ["mkl_rt"]

    libraries = ["lapack", "blas"]
    for lib in libraries:
        if not system_info.get_info(lib):
            raise LinkError("Cannot find {} library".format(lib))
    return libraries


def configuration(parent_package="", top_path=None):
    import pathlib
    from numpy.distutils.misc_util import Configuration

    config = Configuration("sfv", parent_package, top_path)

    fsource = pathlib.Path(__file__).parents[1] / "lib" / "sfv.f90"

    f2pyOnly = [
        "only:",
        "sfv_v",
        "predict_spatial_flux",
        ":",
    ]
    config.add_extension(
        "_sfv",
        sources=[str(fsource), "sfv/_sfv.pyf"],
        libraries=getlibraries(),
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        f2py_options=f2pyOnly,
    )

    return config


if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(**configuration(top_path="").todict())
