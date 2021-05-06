from setuptools import setup

setup(name='milopy',
      version='0.0.99',
      description='python implementation of miloR for differential abundance analysis in single-cell datasets',
      url='https://github.com/emdann/milopy',
      author='Emma Dann',
      author_email='ed6@sanger.ac.uk',
      license='MIT',
      packages=['milopy'],
      install_requires=[
          "pandas",
          "anndata",
          "scanpy",
          "scipy",
          "numpy",
          "matplotlib"
      ],
      zip_safe=False)