# FROM ubuntu:18.04
FROM ensemblorg/ensembl-vep:release_102.0
USER root

WORKDIR /

RUN apt-get update && apt-get install -y vim wget curl git python2.7

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p /miniconda
ENV PATH="/miniconda/bin:$PATH"

# get env file seperately so it doesn't reinstall every time
COPY ./env.yaml /pecgs-somaticwrapper/env.yaml
RUN conda env create --file /pecgs-somaticwrapper/env.yaml 

RUN ln -s /usr/bin/python2.7 /usr/bin/python2

COPY . /pecgs-somaticwrapper

ENV PATH="/usr/bin:/miniconda/envs/somaticwrapper/bin:$PATH"

CMD /bin/bash
