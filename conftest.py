# We never want to do the pybind11 unit tests, so make pytest ignore them.
from pathlib import Path

TO_IGNORE = "pybind11/test"

def pytest_ignore_collect(collection_path, path, config):
    # suppose our condition is some command line argument is passed in
    val = config.getvalue("-k")
    if val == "":
        stringy = collection_path.absolute().as_posix()
        if "cpp_sgwj/pybind11" in stringy:
            return True
    return False

