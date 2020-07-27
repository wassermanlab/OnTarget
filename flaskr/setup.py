import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="OnTargetServer",
    version="0.0.1",
    author="Matthew Mong",
    author_email="mmong@cmmt.ubc.ca",
    description="Lightweight backend for OnTarget",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wassermanlab/OnTarget",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)