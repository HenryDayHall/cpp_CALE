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
            if name.startswith("sgwj_pybind") and name.endswith(".so"):
                return build_dir
    # This file is in the location of the project root
    this_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(build_dir)
    if force_rebuild:
        os.system("make clean")
    os.system("cmake " + this_dir)
    os.system("make -j4")
    return build_dir


def get_module(build_dir):
    # add the build directory to the PYTHONPATH
    sys.path.append(build_dir)
    # import the python bindings
    import sgwj_pybind
    return sgwj_pybind
