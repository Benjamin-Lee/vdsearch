from setuptools import setup, find_packages
import nimporter

setup(
    name="vdsearch",
    version="0.0.1",
    packages=["vdsearch"],
    # This is all the effort required to bundle all Nim modules/libraries
    ext_modules=nimporter.build_nim_extensions(danger=True),
    entry_points={
        "console_scripts": [
            "vdsearch = vdsearch.main:app",
        ],
    },
)
