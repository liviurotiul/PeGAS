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
    version='1.2.3',
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "pegas=pegas.main:main",
            "pegas-lite=pegas.main_lite:main",
        ],
    },
    cmdclass={
        'install': CustomInstallCommand,
    },
    data_files=[
        ('', ['install.sh']),  # Include the install script
    ],
)
