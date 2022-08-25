from setuptools import setup

setup(name='milopy',
      version='0.1.0',
      description='python implementation of miloR for differential abundance analysis in single-cell datasets',
      url='https://github.com/emdann/milopy',
      author='Emma Dann',
      author_email='ed6@sanger.ac.uk',
      license='MIT',
      packages=['milopy'],
      install_requires=[
          "pandas",
          "anndata",
          "scanpy>=1.6.0",
          "mudata>=0.2.0",
          "scipy",
          "numpy",
          "matplotlib",
          "rpy2 >= 3.3.5",
          'leidenalg'
      ],
      zip_safe=False)
