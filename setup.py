from setuptools import find_packages, setup

setup(
    name="vdsearch",
    version="0.0.2",
    packages=find_packages(),
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
        "pycirclize",
        "matplotlib",
    ],
    include_package_data=True,
)
