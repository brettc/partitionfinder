from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import sys

# This file is only used for compiling tigger, so let's just make it do this
# by default. This only works for Mac right now.
sys.argv = ['build_ext.py', 'build_ext', '--inplace']

# CFLAGS
# -Wno-unused-function
extensions = [
    Extension(
        "partfinder/_tiger",
        ["partfinder/_tiger.pyx"],
        extra_compile_args=[
            '-Wno-unused-function',
            '-stdlib=libc++',
            '-std=c++11',
            '-mmacosx-version-min=10.8',
            '-I/usr/local/include',
        ],
        language='c++',
    )
]
setup(
    name='partfinder',
    ext_modules=cythonize(
        extensions,
        aliases={'NUMPY_PATH': numpy.get_include()},
        nthreads=4,
    ),
)
