#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2023/4/11 11:48
# @Author  : zhangchao
# @File    : setup.py
# @Software: PyCharm
# @Email   : zhangchao5@genomics.cn
import setuptools
from wheel.bdist_wheel import bdist_wheel

__version__ = "1.0.0"

class BDistWheel(bdist_wheel):
    def get_tag(self):
        return (self.python_tag, "none", "linux_X86_64")

cmdclass = {
    "bdist_wheel": BDistWheel,
}

requirements = open("requirements.txt").readline()

setuptools.setup(
    name="BatchQC",
    version=__version__,
    author="zhangchao",
    author_email="zhangchao5@genomics.cn",
    url="https://github.com/BGIResearch/batchqc-pipeline.git",
    description="BatchEval",
    python_requires=">=3.8",
    packages=setuptools.find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    cmdclass=cmdclass,
    package_data={'': ["*.html"]},
)


