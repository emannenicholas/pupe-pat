FROM docker.lco.global/docker-miniconda3:4.4.10
MAINTAINER Las Cumbres Observatory <webmaster@lco.global>

RUN conda install -y numpy astropy pytest pytest-runner ipython matplotlib scipy \
        && conda install watchdog emcee corner -c conda-forge \
        && conda clean -y --all

RUN yum -y install epel-release gcc \
        && yum install -y fpack \
        && yum -y clean all

RUN mkdir /home/eng \
        && /usr/sbin/groupadd -g 500 "eng" \
        && /usr/sbin/useradd -g 500 -d /home/eng -M -N -u 500 eng \
        && chown -R eng:eng /home/eng

RUN pip install lcogt-logging Cython \
        && rm -rf ~/.cache/pip

RUN git clone https://github.com/mstamy2/PyPDF2 /usr/src/pypdf2
WORKDIR /usr/src/pypdf2
RUN python setup.py install

RUN git clone https://github.com/kbarbary/sep.git /usr/src/sep
WORKDIR /usr/src/sep
RUN python setup.py install

RUN git clone https://github.com/astropy/astroscrappy.git /usr/src/astroscrappy
WORKDIR /usr/src/astroscrappy
RUN python setup.py install

COPY . /pupepat/src/

WORKDIR /pupepat/src

RUN python setup.py install

ENV HOME /home/eng

WORKDIR /home/eng

USER eng
