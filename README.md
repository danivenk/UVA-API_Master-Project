<img src="https://github.com/danivenk/UVA-API_Master-Project/assets/8749733/ae75f07a-f24a-4a43-b9bb-c2613d3c8c9e" alt="UvA_logo" height="150" />
<img src="https://github.com/danivenk/UVA-API_Master-Project/assets/8749733/f4521512-3629-4472-8106-d311a53e70c7" alt="VU_logo" height="150" />
<img src="https://github.com/danivenk/UVA-API_Master-Project/assets/8749733/12533b55-0573-479d-a5fc-3529247d8570" alt="API_logo" height="150" />

This repository contains the files used for the Master Project at the Anton Pannekoek Institute & University of Amsterdam.
The thesis of this repository and its findings can be found at the [UvA Scriptie Bibliotheek](https://scripties.uba.uva.nl/search?id=record_54076).

# Abstract
Extreme objects in the universe have been found to be in multi-object systems like binaries. When one of these components is a black hole, these binaries can emit X-rays from its accretion disk and the area closer to and around the black hole called the corona. However, due to a still debated mechanism, variations in the X-rays from the accretion disk and the corona will lag behind each other. Recently, a new model has come out to describe these lags. Initially, the model was created in Python, and this project aimed to look into the viability of converting it to a compiled language like C++. The reasons are quicker run time and the ability to run in a standardized fitting program, XSPEC. The latter, however, was not investigated in this thesis. While most of the functions part of this Python lag model could be easily converted with little to no changes, the calculations concerning the illumination fractions could be overhauled entirely by separating the radial dependent and independent calculations. The overhaul caused a speedup of this calculation of around 40-70% depending on the grid size of the geometry used. This speedup is equal to a speedup of several seconds on a function that takes several seconds. To test the converted functions, the C++ version and Python version of the model were fitted to XMM Newton data from GX 339-4 and NICER data from MAXI J1820+070. Both versions of the model fitted the data with similar results and fit quality. In conclusion, the conversion does speed up the model run time but less than one might expect. As such, it remains a question whether it is worth converting the model for use in the standardized fitting program or just implementing the improvements to the Python version.

# Contents
This repository contains the following:
- All the images and figures used for the thesis as well as the thesis versions.
- A test for the creation of a pybind module as well as trying out XSPEC in `test/` & `pymodule/`
- The flowcharts of all the functions in the model in `flowcharts/`.
- The original model as well as some of the fitting notebooks in `mono_lags/`
- The converted model to C++ in `C++/`

## `mono_lags/`
In this directory all the files related to the original model and fitting this model that were acquired. This also includes two datasets from GX 339-4 ([Uttley et. al. 2011](https://ui.adsabs.harvard.edu/abs/2011MNRAS.414L..60U/)) & MAXI J1820+070 ([Wang et. al. 2021](https://ui.adsabs.harvard.edu/abs/2021ApJ...910L...3W); [Bollemeijer et. al. 2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.528..558B)) to check if the converted model and the original model would get similar results (proving the conversion worked). It uses modules like [Numpy](https://numpy.org/), [Scipy](https://scipy.org/).

## `C++/`
In this C++ directory all the files related to the converted model as well as the C++ side of testing. It converted the model that is located in `mono_lags/`. For sharing values between the versions it uses [nlohmann/json](https://github.com/nlohmann/json) and for the fast fourier transform [mreineck/pocketfft](https://github.com/mreineck/pocketfft).
