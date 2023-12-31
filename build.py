import os
import sys


def build(build_dir = "./build", force_rebuild = False):
    # get the absolute path of the build directory
    build_dir = os.path.abspath(build_dir)
    if not os.path.exists(build_dir):
        os.mkdir(build_dir)
    elif not force_rebuild:
        # check if the python bindings already exist
        for name in os.listdir(build_dir):
            if name.startswith("CALE_pybind") and name.endswith(".so"):
                return build_dir
    # This file is in the location of the project root
    this_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(build_dir)
    if force_rebuild:
        os.system("make clean")
    os.system("cmake " + this_dir)
    os.system("make -j4")
    return build_dir

CALE_pybind = None


def get_module(build_dir):
    global CALE_pybind
    if CALE_pybind is not None:
        return CALE_pybind
    # add the build directory to the PYTHONPATH
    sys.path.append(build_dir)
    # import the python bindings
    import CALE_pybind
    return CALE_pybind
