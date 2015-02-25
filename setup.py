from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import sys

# HACK for now
sys.argv = ['setup.py', 'build_ext', '--inplace']

# CFLAGS
# -Wno-unused-function
extensions = [
    Extension(
        "partfinder/_tiger", 
        ["partfinder/_tiger.pyx"], # "bricolage/scheme_cooperative.cpp"],
        extra_compile_args = [
            # '-ffast-math',
            '-Wno-unused-function', 
            # '-stdlib=libc++',
            # '-std=c++11', 
            '-mmacosx-version-min=10.8',
        ],
        # depends = ['bricolage/pubsub2_c.h'],
        language = 'c++',
        # include_path = [numpy.get_include()],
        # include_dirs = ['./bricolage'],
        # libraries = ['bricolage'],
        # library_dirs = ['./bricolage'],
    )
]
setup(
    name='partfinder',
    ext_modules=cythonize(
        extensions,
        aliases={ 'NUMPY_PATH': numpy.get_include() },
        nthreads=4,
    ),
)
