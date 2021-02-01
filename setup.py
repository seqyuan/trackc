import os, re
from setuptools import setup, find_packages
 
with open("README.md", "r") as fh:
    long_description = fh.read()

def get_version():
    try:
        f = open(os.path.join("trackc", "_version.py"))
    except EnvironmentError:
        return None
    for line in f.readlines():
        mo = re.match("__version__ = '([^']+)'", line)
        if mo:
            ver = mo.group(1)
            return ver
    return None

setup(
    name="trackC",
    version=get_version(),
    author="ahworld",
    author_email="seqyuan@gmail.com",
    description="Base functions to plot Xomics data tracks",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/seqyuan/trackC",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)




