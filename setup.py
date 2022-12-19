"""
OnTarget: in silico design of MiniPromoters for targeted delivery of expression
"""

import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="ontarget",
    version="22.12.1",
    author="Oriol Fornes, Tamar Av-Shalom",
    author_email="oriol.fornes@gmail.com, avshalom.tamar0@gmail.com",
    description="OnTarget: in silico design of MiniPromoters for targeted delivery of expression",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wassermanlab/OnTarget",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
