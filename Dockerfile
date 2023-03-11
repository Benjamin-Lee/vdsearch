
FROM nimlang/nim:1.6.10

# install dependencies
RUN apt update && apt install -y \
    python3 \
    python3-pip \
    mmseqs2 \
    infernal \
    # for RNAmotif
    flex \
    bison \
    # for ViennaRNA
    libgsl23

RUN pip install nimporter==1.1.0
RUN nimble refresh && nimble install cello && nimble install nimpy

# install seqkit
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.3.1/seqkit_linux_amd64.tar.gz -O ~/seqkit.tar.gz && \
    tar xf ~/seqkit.tar.gz && \
    cp seqkit /usr/local/bin/

# install ViennaRNA
RUN wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_20_04/viennarna_2.5.1-1_amd64.deb && \
    dpkg -i viennarna_2.5.1-1_amd64.deb && \
    apt --fix-broken install -y

# install RNAmotif
RUN git clone https://github.com/dacase/rnamotif && \
    cd rnamotif && \
    make && \
    cp src/rnamotif /usr/bin/ && \
    cp src/rmprune /usr/bin/ && \
    cp src/rmfmt /usr/bin/

WORKDIR /vdsearch
COPY ./vdsearch ./vdsearch
COPY ./setup.py .
COPY ./MANIFEST.in .
RUN pip install .

# download viroiddb and CM database since the image has its own filesystem
RUN vdsearch download-viroiddb && \
    vdsearch download-cms
 
WORKDIR /data
ENTRYPOINT ["vdsearch"]