from setuptools import setup, find_packages

# Function to read the requirements.txt file
def read_requirements():
    with open('requirements.txt') as req_file:
        return req_file.read().splitlines()

setup(
    name="pistacker",
    author="Eric D. Sakkas",
    author_email="esakkas@wesleyan.edu",
    description="Analyzes Pi-Stacking in Molecular Dynamics Trajectories",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/esakkas24/stacker",
    version="1.1.4",
    packages=find_packages(),
    python_requires=">=3.10",
    entry_points={
        'console_scripts': [
            'stacker=stacker.__init__:run_python_command', 
        ],
    },
    install_requires=read_requirements(), 
)