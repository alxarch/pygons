from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

pygons_src = ['src/mathutils.c',\
			  'src/bounds.c',\
			  'src/polygon.c',\
			  'src/split.c',\
			  'src/boolean.c',\
			  'src/pygons.pyx']

setup(
    name = 'pyGons',
    author = 'Alexandros Sigalas',
    author_email = 'alxarch@gmail.com',
    version = '0.1',
    description = 'A polygon boolean operations module for Python implemented in C.',
    license = 'GNU General Public License (GPL)',
    packages = ['pygons'],
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension('pygons', pygons_src, libraries=["m"])]
)
