import nimporter
from setuptools import find_packages, setup

setup(
    name="vdsearch",
    version="0.0.1",
    packages=find_packages(),
    # This is all the effort required to bundle all Nim modules/libraries
    ext_modules=nimporter.build_nim_extensions(danger=True),
    entry_points={
        "console_scripts": [
            "vdsearch = vdsearch.main:app",
        ],
    },
    install_requires=[
        "click-didyoumean",
        "pandas",
        "rich",
        "rich_click",
        "shellingham",
        "typer<0.6",
        "scikit-bio",
        "igraph",
    ],
    include_package_data=True,
)
