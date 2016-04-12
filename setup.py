from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='quaternion',
      version='0.1',
      description='A package to handle quaternions',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
      ],
      keywords='quaternion transformations',
      url='https://github.com/mjsobrep/python_quaternions',
      author='Michael Sobrepera',
      author_email='mjsobrep@live.com',
      license='Copyright (c) 2016 GTRC',
      packages=['quaternion'],
      install_requires=[],
      include_package_data=True,
      zip_safe=False,
      test_suite="quaternion.tests"
      )
