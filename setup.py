import setuptools
import glob
import os
#from setuptools import setup, find_packages


def read_requirements(fname):
    with open(fname, 'r', encoding='utf-8') as file:
        return [line.rstrip() for line in file]


setuptools.setup(
    name='trackc',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    packages=setuptools.find_packages(where='src'),
    package_dir={'': 'src'},
    py_modules=[os.path.splitext(os.path.basename(path))[0] for path in glob.glob('src/*.py')],
    include_package_data=True,
    install_requires=read_requirements('requirements.txt'),
    author="Zan Yuan",
    author_email="yfinddream@gmail.com",
    #packages=find_packages(),
    scripts=['bin/gtf2bed4trackc']
    #include_package_data=True,
    #package_dir={'trackc': 'src/trackc'},
    #package_data={'trackc': ['qc_template.html']},
    url='http://trackc.readthedocs.io',
            
    description="Track view of chromosome conformation and multi-omics data",
    long_description=open('README.rst').read(),
    #url="https://github.com/seqyuan/trackc",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
