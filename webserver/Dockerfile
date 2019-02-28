FROM ubuntu:latest

MAINTAINER Zealseeker "yanyanghong@163.com"


# Get conda, supervisor
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 WEBENV=docker
RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git mercurial subversion libxrender-dev supervisor build-essential \
	&& apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN echo "export PATH=/opt/conda/bin:$PATH" > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm ~/anaconda.sh

# Set and create environment
ENV PATH /opt/conda/bin:$PATH

RUN conda install -y -c rdkit rdkit flask && \
    conda install -y -c conda-forge uwsgi

COPY . /var/www/chemviewer/
CMD /usr/bin/supervisord -c /var/www/chemviewer/webserver/supervisord.conf
