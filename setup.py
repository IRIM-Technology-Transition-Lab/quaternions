from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='dual_quaternions',
      version='0.1',
      description='A package to handle dual quaternions',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Manufacturing',
        'Intended Audience :: Science/Research',
        'License :: Other/Proprietary License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
      ],
      keywords='quaternions robotics transformations',
      url='https://github.com/mjsobrep/python_dual_quaternions',
      author='Michael Sobrepera',
      author_email='mjsobrep@live.com',
      license='Copyright (c) 2016 Michael Sobrepera',
      packages=['ur_cb2', 'ur_cb2.receive', 'ur_cb2.send'],
      install_requires=[
      ],
      include_package_data=True,
      entry_points={
        'console_scripts': ['cb2-listen=ur_cb2.receive.cb2_receive_example:'
                            'main',
                            'cb2-record=ur_cb2.receive.cb2_store_points:main',
                            'cb2-play=ur_cb2.cb2_move_to_points:main']
      },
      zip_safe=False)
