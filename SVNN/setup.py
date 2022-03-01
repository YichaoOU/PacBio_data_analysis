from setuptools import setup, find_packages

setup(name='svnn',
      version='0.1.0',
      description='SVNN: Structural Variation Caller',
      author='Shaya Akbarinejad',
      author_email='akbarinejadshaya@gmail.com',
      zip_safe=False,
      install_requires=['pysam', 'numpy', 'scipy', 'matplotlib'],
      scripts=['source/svnn']
     )