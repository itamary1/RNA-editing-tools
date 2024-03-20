############################## Dockerfile ###################################

FROM centos:centos7.5.1804

###################################### METADATA ##############################

################## BEGIN INSTALLATION ######################
#
# Change to root
#
USER root
RUN echo "update OS"
RUN yum -y update && yum clean all
RUN echo "install Libs"
RUN yum install gcc \
    openssl-devel \
    bzip2-devel \
    libffi-devel \
    zlib-devel \
    dnf \
    sed \
    curl \
    fuse \
    git \
    wget \
    zip \
    pkg-config \
    bzip2 \
    gzip \
    gawk \
    grep \
    git \
    sudo \
    ca-certificates \
    -y && yum clean all
RUN dnf update -y
RUN dnf install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm -y
RUN sudo dnf makecache --refresh
RUN yum group install -y "Development Tools" && yum clean all



#
# INSTALL bioconda and python
#
RUN curl -sSL https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local/ \
    && rm -rf /tmp/miniconda.sh 
RUN conda install -y python=2.7 \
    && conda update conda \
    && conda clean --all --yes \
    && rpm -e --nodeps curl bzip2 \
    && yum clean all
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

#RUN conda install --yes --channel=cyclus java-jdk=8.0.92

#
# Add biodocker user with password biodocker
#
RUN useradd --create-home --shell /bin/bash --user-group --uid 1000 biodocker && \
    echo `echo "biodocker\nbiodocker\n" | passwd biodocker` && sudo usermod -aG wheel biodocker




#
# Create $HOME/bin folder
#
RUN mkdir -p /home/biodocker/bin && chmod -R 777 /home/biodocker

RUN cd /home/biodocker &&\
    git clone https://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@github.com/shalomhillelroth/GGPS




RUN chmod -R 777 /home/biodocker/
#
# Add $HOME/bin to path
#
ENV PATH=$PATH:/home/biodocker/bin
ENV HOME=/home/biodocker

#
# INSTALL SAMTOOLS and BEDTOOLS
#
#RUN conda install -y bedtools
RUN conda install -c bioconda bedtools -y

RUN dnf install samtools -y

#
#INSTALL bamUtils
#
RUN mkdir bumUtilInstall && cd bumUtilInstall && \
    git clone https://github.com/statgen/bamUtil && cd bamUtil/ && make cloneLib && \
    make && make install INSTALLDIR=/home/biodocker/bin && cd .. && rm -rf ./bumUtilInstall


RUN dnf install java-11-openjdk -y

#
# editingIndexUtils
#
ADD ./EditingIndexJavaUtils.jar /home/biodocker/bin
ADD ./PileupKnown.conf  /home/biodocker/bin
ADD ./PileupKnown.SoftIntersect.conf  /home/biodocker/bin



# Change user
#
# USER biodocker

#
# Add $HOME/bin to path
#
ENV PATH=$PATH:/home/biodocker/bin
ENV HOME=/home/biodocker




#
# Share default volumes
#
VOLUME ["/data", "/config"]



CMD ["/usr/bin/bash"]

##################### INSTALLATION END #####################

################################### TRY SOLVE NEXTFLOW USER PROBLEM

# # Create a new user and group inside the container
# RUN groupadd -r myusergroup && useradd -r -g myusergroup myuser

# # Set the working directory and ownership
# WORKDIR /app
# RUN chown -R myuser:myusergroup /app

# # Set the user to be used when the container starts
# USER myuser
