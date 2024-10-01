import os
from setuptools import setup, find_packages
from setuptools.command.install import install
import subprocess

class CustomInstallCommand(install):
    def run(self):
        install.run(self)
        subprocess.call(['bash', 'install.sh'])

setup(
    name='pegas',
    version='0.1.14',
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    cmdclass={
        'install': CustomInstallCommand,
    },
    data_files=[
        ('', ['install.sh']),  # Include the install script
    ],
    entry_points={
        'console_scripts': [
            'pegas=pegas.main:main',
        ],
    },
)
