FROM ubuntu

ARG TRAVIS=true
ARG TRAVIS_JOB_ID
ARG PYTHON_VERSION=3

# Copy latest source code into the image
WORKDIR /test
WORKDIR pyensembl
ADD . .

# System setup
RUN apt-get update
# wget, to be able to download Conda script
RUN apt-get install -y wget
# bzip2, so that Conda can unbundle stuff
RUN apt-get install -y bzip2
# because http://stackoverflow.com/a/25423366
RUN mv /bin/sh /tmp/sh && ln -s /bin/bash /bin/sh

### Conda Setup ###
WORKDIR ../programs
RUN wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -f -p ./miniconda

## pyensemble setup and tests
WORKDIR ../pyensembl
RUN export PATH="`pwd`/../programs/miniconda/bin:$PATH" && \
    # Create our Conda environment and activate it
    conda create -y -n pyensembl-test python=${PYTHON_VERSION} && \
    . activate pyensembl-test && \

    # Install some packages up-front
    conda install -y numpy scipy nose pandas matplotlib biopython && \

    # and install requirements
    pip install coveralls && \
    pip install -r requirements.txt && \

    # and install pyensembl
    pip install . && \

    # Cache all the genomes!
    pyensembl install --release 54 --species human && \
    pyensembl install --release 75 --species human && \
    pyensembl install --release 77 --species human && \
    pyensembl install --release 83 --species human && \
    pyensembl install --release 84 --species human && \
    pyensembl install --release 67 --species mouse && \
    pyensembl install --release 84 --species mouse && \

    # test it
    ./lint.sh && \
    nosetests test --verbose --with-coverage --cover-package=pyensembl && \
    if [ "$TRAVIS" == "true" ]; then coveralls; fi

## Clean up to reduce the image size
# Empty test folders
RUN rm -rf /test
# Remove additional software
RUN apt-get remove -y wget bzip2
RUN rm -rf ~/.cache/pip
# Empty optional/downloaded genome files
RUN export CACHE=~/.cache/pyensembl && \
    find ~/.cache/pyensembl -name '*.csv' -exec sh -c > {}" \; && \
    find ~/.cache/pyensembl -name '*.fa.gz' -exec sh -c > {}" \;
# Put the standard sh back in
RUN mv -f /tmp/sh /bin/sh

## All done!
RUN echo "Latest version of pyensembl successfully built and tested on Python ${PYTHON_VERSION}.x"