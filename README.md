# Spectral Graph Wavelet Jet clustering

To compile;

- Ensure you have pulled the pybind11 submodule. If you haven't then call `git submodule update --init --recursive`, to get pybind11.
- Make a build dir somewhere outside the repo, and enter the build dir.
- From the build dir, call `cmake ../path/to/base/of/repo/ && make`.
- Check the C++ worked with `./CALE_toy`.
- Check the python bindings work with `ipython3 -c "import CALE_pybind; CALE_pybind.PxPyPz(6., 3., 1., 2.)"`.
- ???
- Profit.

