import os
import re
from distutils.core import setup


def version():
    setup_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(setup_dir, 'src', 'VERSION')) as fh:
        return re.search(r'software_version=(.+)\n', fh.read()).group(1)


def add_py_scripts(root):
    out = list()
    for subdir, dirs, files in os.walk(root):
        for cur_f in files:
            cur_path = os.path.join(subdir, cur_f)
            if cur_f != '__init__.py':
                out.append(cur_path)
            print os.path.join(subdir, cur_f)
    return out


py_scripts = ['bin/gtdb']
py_scripts.extend(add_py_scripts('scripts_user'))
py_scripts.extend(add_py_scripts('scripts_dev'))

setup(
    name='gtdb',
    version=version(),
    author='Pierre-Alain Chaumeil and Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['gtdb'],
    package_dir= {'gtdb': 'src/gtdb'},
    scripts=py_scripts,
    package_data={'gtdb': ['VERSION',
                           'MANIFEST.in']},
    data_files=[('', ['src/VERSION', 'README.md'])],
    url='https://github.com/Ecogenomics/GTDBNCBI',
    license='GPL3',
    description='The GTDB provides the software infrastructure for working with a large collections of genomic resources.',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.8.0",
        'biolib < 0.1.0',
        'requests',
        'unidecode',
        'matplotlib',
        'python-dateutil'
    ]
)
