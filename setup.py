# to compile cython sources, call:
# python setup.py install
#
# Compilation can otherwise be obtained via the following commands:
# cython -a graph.pyx
# gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.6 -o graph.so graph.c
# cython -a fastmath.pyx
# gcc -shared -pthread -fPIC -fwrapv -O3 -Wall -fno-strict-aliasing -I/usr/include/python2.6 -o fastmath.so fastmath.c

import os, stat, sys
import shutil
import numpy as np
from distutils.core import setup
from distutils.command.build_ext import build_ext
from Cython.Build import cythonize

class InstallCommand(build_ext):

    def run(self):

        build_ext.run(self)

        try:
            for root, dirnames, filenames in os.walk("build"):
                print dirnames
                for filename in filenames:
                    extension = filename.split(".")[-1]
                    if extension in ["pyd", "so"]:
                        os.rename(os.path.join(root, filename), filename)

        except Exception as ex:
            print("files already exist, skipping...")

        shutil.rmtree("build")

# Cythonising necessary scripts function
def cythonise(folder):
    """
    Cythonise pyx files 

    :param folder: Location of pyx files
    """

    home = os.getcwd()
    os.chdir(folder)

    # small hack to get around a problem in older cython versions, i.e.
    # an infinite dependencies loop when __init__.py file is in the same folder as pyx
    if os.path.exists("__init__.py"):
        os.rename("__init__.py", "tmp")

    setup(
        include_dirs=[np.get_include()],
        ext_modules=cythonize(
            "*.pyx",
            include_path=[np.get_include()],
            compiler_directives={'boundscheck': False, 'wraparound': False, 'embedsignature': True}),
        cmdclass={'install': InstallCommand}
    )

    # continuation of the small hack
    os.rename("tmp", "__init__.py")

    os.chdir(home)


# Add auto_scripts to path
home = os.getcwd()
sys.path.append(os.getcwd())
os.chdir("auto_scripts")

for filename in os.listdir(os.getcwd()):
    if filename.endswith(".sh") or filename.endswith(".py"):
        os.chmod(filename, stat.S_IRWXU)
    else:
        continue

os.chdir(home)

# Cythonise files in methods
cythonise('methods/')
