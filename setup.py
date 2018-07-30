# Copyright (c) 2014-2018. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function
import os
import re

from setuptools import setup

package_name = "pyensembl"
current_directory = os.path.dirname(__file__)
readme_filename = 'README.md'
readme_path = os.path.join(current_directory, readme_filename)
github_url = "https://github.com/openvax/%s" % package_name

try:
    with open(readme_path, 'r') as f:
        readme_markdown = f.read()
except IOError as e:
    print(e)
    print("Failed to open %s" % readme_path)
    readme_markdown = ""

try:
    import pypandoc
    readme_restructured = pypandoc.convert(readme_markdown, to='rst', format='md')
except ImportError as e:
    readme_restructured = readme_markdown
    print(e)
    print("Failed to convert %s to reStructuredText" % readme_filename)
    pass

with open('%s/__init__.py' % package_name, 'r') as f:
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        f.read(),
        re.MULTILINE).group(1)

if not version:
    raise RuntimeError('Cannot find version information')

if __name__ == '__main__':
    setup(
        name=package_name,
        version=version,
        description="Python interface to ensembl reference genome metadata",
        author="Alex Rubinsteyn",
        author_email="alex.rubinsteyn@mssm.edu",
        url=github_url,
        license="http://www.apache.org/licenses/LICENSE-2.0.html",
        entry_points={
            'console_scripts': [
                'pyensembl = %s.shell:run' % package_name
            ],
        },
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache Software License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        install_requires=[
            "typechecks>=0.0.2",
            "numpy>=1.7",
            "pandas>=0.15",
            "datacache>=1.1.4",
            "memoized-property>=1.0.2",
            "six>=1.9.0",
            "gtfparse>=1.1.0",
            "serializable",
            "tinytimer",
        ],
        long_description=readme_restructured,
        packages=[package_name],
        package_data={package_name: ['logging.conf']},
    )
