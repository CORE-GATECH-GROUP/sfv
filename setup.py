from numpy.distutils.core import setup


def configuration(parent_package="", top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package=parent_package, top_path=top_path)
    config.add_subpackage("sfv")

    return config


SETUP_OPTS = dict(
    version="0.3.2",
    name="sfv",
    author="Andrew Johnson",
    author_email="ajohnson400@gatech.edu",
    url="https://pwp.gatech.edu/core",
    install_requires=["wheel", "setuptools", "numpy>=1.7", "scipy"],
)

SETUP_OPTS["configuration"] = configuration

setup(**SETUP_OPTS)
