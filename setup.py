import ez_setup
import sys
ez_setup.use_setuptools()
from setuptools import setup

# from mpld3
def get_version(path):
    """Get the version info from package without importing it"""
    import ast

    with open(path) as init_file:
        module = ast.parse(init_file.read())

    version = (ast.literal_eval(node.value) for node in ast.walk(module)
               if isinstance(node, ast.Assign)
               and node.targets[0].id == "__version__")
    try:
        return next(version)
    except StopIteration:
        raise ValueError("version could not be located")


setup(name='cigar',
      version=get_version("cigar.py"),
      description="manipulate SAM cigar strings",
      py_modules=['cigar'],
      author="Brent Pedersen",
      author_email="bpederse@gmail.com",
      url="https://github.com/brentp/cigar",
      license="MIT",
      long_description=open('README.md').read(),
      classifiers=[
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'Programming Language :: Python :: 2',
      'Programming Language :: Python :: 3'
      ],
)
