from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext, \
    ParallelCompile, naive_recompile

ext_modules = [
    Pybind11Extension(
        "example",
        sorted(glob("*.cpp")),
    ),
]

ParallelCompile("NPY_NUM_BUILD_JOBS", needs_recompile=naive_recompile).install()
setup(name="my_ex", version="0.0.1", cmdclass={"build_ext": build_ext}, ext_modules=ext_modules)