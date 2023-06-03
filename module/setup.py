from numpy.distutils.core import setup, Extension

module_f = Extension(name="fit_ring",sources=["src/fit_ring.f90"]) 
# name: where the .so file are created

setup(
    name="fit_ring", # name of the module
    version="0.1",
    ext_modules=[module_f],
    description="Fit observations with axisymmetric ring model",
    author="Kiyoaki Doi",
)

