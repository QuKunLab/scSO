from setuptools import setup
from setuptools import Extension

SNMF_inner_module = Extension(name='SNMF_inner',  # 模块名称
                              sources=['SNMF_inner.cpp'],    # 源码
                              include_dirs=[r'..\eigen-3.3.8',
                                            r'.\env\Scripts',     # 依赖的第三方库的头文件
                                            r'.\env\lib\python3.8\site-packages\pybind11']
                              )

setup(ext_modules=[SNMF_inner_module])
