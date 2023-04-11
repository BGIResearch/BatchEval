#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2023/4/11 11:48
# @Author  : zhangchao
# @File    : setup.py
# @Software: PyCharm
# @Email   : zhangchao5@genomics.cn
from setuptools import setup, find_packages

__version__ = "1.0.0"

requirements = open("requirements.txt").readline()

setup(
    name="BatchQC",
    version=__version__,
    author="zhangchao",
    author_email="zhangchao5@genomics.cn",
    url="https://github.com/BGIResearch/batchqc-pipeline.git",
    description="BatchQC Pipeline",
    python_requires=">=3.8.8",
    packages=find_packages(),
    install_requires=requirements
)


