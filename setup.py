from setuptools import setup

with open("requirements.txt") as f:
    required_packages = f.readlines()

required_packages = [r for r in required_packages if "-e" not in r]

setup(
    name="lc_correction",
    version="1.0.0",
    description='Scripts for ALeRCE light curve correction',
    author="ALeRCE Team",
    author_email='contact@alerce.online',
    packages=['lc_correction'],
    install_requires=required_packages,
    build_requires=required_packages
)