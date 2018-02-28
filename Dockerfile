FROM docker.lcogt.net/miniconda3:4.3.21
MAINTAINER Las Cumbres Observatory <webmaster@lco.global>

RUN conda install -y numpy astropy pytest pytest-runner ipython matplotlib scipy \
        && conda install sep -c openastronomy \
        && conda install watchdog -c conda-forge \
        && conda clean -y --all

RUN yum -y install epel-release \
        && yum install -y fpack \
        && yum -y clean all

RUN mkdir /home/eng \
        && /usr/sbin/groupadd -g 500 "eng" \
        && /usr/sbin/useradd -g 500 -d /home/eng -M -N -u 500 eng \
        && chown -R eng:eng /home/eng

RUN pip install lcogt-logging \
        && rm -rf ~/.cache/pip

COPY . /pupepat/src/

WORKDIR /pupepat/src

RUN python setup.py install

ENV HOME /home/eng

WORKDIR /home/eng
