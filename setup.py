from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="pybioanalyzer",
    version="0.1",
    description="A Python parser for Bioanalyzer files",
    author="Shunsuke Sumi, Evgeniia Edeleva, Georg Urtel",
    author_email="shunsuke.sumi@cira.kyoto-u.ac.jp",
    url="https://github.com/Shunsuke-1994/PyBioAnalyzer",
    license="MIT",
    packages=["pybioanalyzer"],
    python_requires=">=3.7.0",
    install_requires=[
        "matplotlib",
        "scipy>=1.2.1",
        "pandas"
    ],  # external packages as dependencies
    long_description="README.md",
    long_description_content_type="text/markdown",
)
