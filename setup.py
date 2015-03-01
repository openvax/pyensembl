# Copyright (c) 2014. Mount Sinai School of Medicine
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

import os


readme_filename = os.path.join(os.path.dirname(__file__), 'README.md')
with open(readme_filename, 'r') as f:
    readme = f.read()


try:
    import pypandoc
    readme = pypandoc.convert(readme, to='rst', format='md')
except:
    print("Conversion of long_description from MD to reStructuredText failed...")
    pass


requires_list = []
try:
    from pip.req import parse_requirements
    parsed_reqs = parse_requirements('requirements.txt')
    requires_list = [str(req.req) for req in parsed_reqs]
except ImportError:
    print("Please install pip before proceeding")


from setuptools import setup

if __name__ == '__main__':
    setup(
        name='pyensembl',
        version="0.5.6",
        description="Python interface to ensembl reference genome metadata",
      	author="Alex Rubinsteyn",
        author_email="alex {dot} rubinsteyn {at} mssm {dot} edu",
      	url="https://github.com/hammerlab/pyensembl",
        license="http://www.apache.org/licenses/LICENSE-2.0.html",
        entry_points={
          'console_scripts': [
              'pyensembl = pyensembl.shell:run'
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
        install_requires=requires_list,
        long_description=readme,
        packages=['pyensembl'],
    )
