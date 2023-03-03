FROM cern/cc7-base:latest

USER root
ENV WKDIR=/root
WORKDIR ${WKDIR}


##########
# Basics #
##########

RUN yum -y install which wget tar devtoolset-8-gcc-c++ make cmake3 git mesa-libGL-devel libXmu-devel expat-devel && \
    yum group install -y "X Window System" && \
    yum clean all
RUN ln -s -f /usr/bin/cmake3 /usr/bin/cmake
ENV DEVTOOLS_BINDIR=/opt/rh/devtoolset-8/root/usr/bin
RUN ln -s -f ${DEVTOOLS_BINDIR}/g++ /usr/bin/g++
RUN ln -s -f ${DEVTOOLS_BINDIR}/gcc /usr/bin/gcc
RUN ln -s -f ${DEVTOOLS_BINDIR}/ld /usr/bin/ld
ENV CC=/opt/rh/devtoolset-8/root/usr/bin/gcc
ENV CXX=/opt/rh/devtoolset-8/root/usr/bin/g++


#########
# HEPMC #
#########

ARG HEPMC_VERSION=2.06.11

RUN wget hepmc.web.cern.ch/hepmc/releases/hepmc${HEPMC_VERSION}.tgz && \
    tar -xzf hepmc${HEPMC_VERSION}.tgz
ENV HEPMC_SRC_DIR=${WKDIR}/HepMC-${HEPMC_VERSION}
ENV HEPMC_DIR=${WKDIR}/hepmc
RUN mkdir ${HEPMC_DIR}
RUN cd ${HEPMC_DIR} && \
    cmake ${HEPMC_SRC_DIR} -Dmomentum:STRING=MEV -Dlength:STRING=MM && \
    make -j4
RUN mkdir ${HEPMC_DIR}/include && \
    cp -r ${HEPMC_SRC_DIR}/HepMC ${HEPMC_DIR}/include/


###########
# FASTJET #
###########

ARG FASTJET_VERSION=3.4.0

RUN cd ${WKDIR}
ENV FASTJET_HOME=${WKDIR}/fastjet
RUN mkdir ${FASTJET_HOME}
RUN wget http://fastjet.fr/repo/fastjet-${FASTJET_VERSION}.tar.gz && \
    tar -xzf fastjet-${FASTJET_VERSION}.tar.gz
RUN cd fastjet-${FASTJET_VERSION} && \
    ./configure --prefix=$FASTJET_HOME && \
    make -j4 && \
    make install


###########
# JSONCPP #
###########

RUN yum -y install jsoncpp-devel
RUN cp -r /usr/include/jsoncpp/json /usr/include


###########
# PYTHIA8 #
###########

ARG PYTHIA_VERSION=8306

RUN cd ${WKDIR}
RUN wget https://pythia.org/download/pythia83/pythia${PYTHIA_VERSION}.tgz && \
    tar -xzf pythia${PYTHIA_VERSION}.tgz
ENV PYTHIA8_HOME=${WKDIR}/pythia${PYTHIA_VERSION}
RUN cd ${PYTHIA8_HOME} && \
    make -j4


########
# ROOT #
########

RUN yum -y install root && \
    yum clean all


############
# Clean Up #
############

RUN cd ${WKDIR}
RUN rm -f *tgz *tar.gz


###########
# Geant 4 #
###########

RUN cd ${WKDIR}
RUN wget http://cern.ch/geant4-data/releases/lib4.10.7.p03/Linux-g++8.3.0-CC7.tar.gz
RUN tar -xzf Linux-g++8.3.0-CC7.tar.gz

ENV Geant4_DIR=${WKDIR}/Geant4-10.7.3-Linux
ENV Geant4_InstallBaseDir=${Geant4_DIR}/share/Geant4-10.7.3
ENV Geant4_DataDir=${Geant4_InstallBaseDir}/data
RUN mkdir ${Geant4_DataDir}

RUN wget http://cern.ch/geant4-data/releases/geant4.10.07.p03.tar.gz
RUN tar -xzf geant4.10.07.p03.tar.gz
RUN sed -i 's/GEANT4_USE_OPENGL_X11 "Build Geant4 OpenGL driver with X11 support" OFF/GEANT4_USE_OPENGL_X11 "Build Geant4 OpenGL driver with X11 support" ON/' ${WKDIR}/geant4.10.07.p03/cmake/Modules/G4InterfaceOptions.cmake
ENV G4_BuildDir=${WKDIR}/g4_build
RUN mkdir ${G4_BuildDir} && \
    cd ${G4_BuildDir} && \
    cmake ${WKDIR}/geant4.10.07.p03 && \
    make -j4
RUN rm -rf ${Geant4_DIR}/lib64/*.so*
RUN cp -r ${G4_BuildDir}/BuildProducts/lib64/* ${Geant4_DIR}/lib64/

RUN cd ${Geant4_DataDir} && \
    wget http://cern.ch/geant4-data/datasets/G4NDL.4.6.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4EMLOW.7.13.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4PII.1.3.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4PhotonEvaporation.5.7.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4RadioactiveDecay.5.6.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4RealSurface.2.2.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4SAIDDATA.2.0.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4ABLA.3.1.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4INCL.1.0.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4ENSDFSTATE.2.3.tar.gz && \
#    wget http://cern.ch/geant4-data/datasets/G4PARTICLEXS3.1.1.tar.gz && \
    tar -xzf G4NDL.4.6.tar.gz && \
    tar -xzf G4EMLOW.7.13.tar.gz && \
    tar -xzf G4PII.1.3.tar.gz && \
    tar -xzf G4PhotonEvaporation.5.7.tar.gz && \
    tar -xzf G4RadioactiveDecay.5.6.tar.gz && \
    tar -xzf G4RealSurface.2.2.tar.gz && \
    tar -xzf G4SAIDDATA.2.0.tar.gz && \
    tar -xzf G4ABLA.3.1.tar.gz && \
    tar -xzf G4INCL.1.0.tar.gz && \
#    tar -xzf G4PARTICLEXS3.1.1.tar.gz && \
    tar -xzf G4ENSDFSTATE.2.3.tar.gz
COPY G4PARTICLEXS3.1.1.tar.gz ${Geant4_DataDir}/G4PARTICLEXS3.1.1.tar.gz
RUN cd ${Geant4_DataDir} && tar -xzf G4PARTICLEXS3.1.1.tar.gz
RUN rm -f *tar.gz

ENV G4ABLADATA=${Geant4_DataDir}/G4ABLA3.1
ENV G4INCLDATA=${Geant4_DataDir}/G4INCL1.0
ENV G4LEDATA=${Geant4_DataDir}/G4EMLOW7.13
ENV G4PARTICLEXSDATA=${Geant4_DataDir}/G4PARTICLEXS3.1.1
ENV G4REALSURFACEDATA=${Geant4_DataDir}/RealSurface2.2
ENV G4LEVELGAMMADATA=${Geant4_DataDir}/PhotonEvaporation5.7
ENV G4PIIDATA=${Geant4_DataDir}/G4PII1.3
ENV G4SAIDXSDATA=${Geant4_DataDir}/G4SAIDDATA2.0
ENV G4ENSDFSTATEDATA=${Geant4_DataDir}/G4ENSDFSTATE2.3
ENV G4NEUTRONHPDATA=${Geant4_DataDir}/G4NDL4.6
ENV G4RADIOACTIVEDATA=${Geant4_DataDir}/RadioactiveDecay5.6
ENV G4WORKDIR=${Geant4_DataDir}/root/geant4_workdir

ENV G4LIB_BUILD_SHARED=1
ENV G4UI_USE_TCSH=1
ENV G4LIB_USE_ZLIB=1
ENV G4VIS_USE_OPENGLX=1
ENV G4INCLUDE=${Geant4_DIR}/include/Geant4
ENV G4INSTALL=${Geant4_InstallBaseDir}/geant4make
ENV G4LIB=${Geant4_DIR}/lib64/Geant4-10.7.3
ENV G4SYSTEM=Linux-g++

#######
# SCD #
#######

ENV SCD_BASE_DIR=${WKDIR}/SCD
ENV SCD_SRC_DIR=${SCD_BASE_DIR}/SCD
ENV SCD_BUILD_DIR=${SCD_BASE_DIR}/build
RUN mkdir ${SCD_BASE_DIR} ${SCD_BUILD_DIR}
COPY . ${SCD_BASE_DIR}
# RUN cd ${SCD_BASE_DIR} && git checkout -b SuperRes origin/SuperRes
RUN cd ${SCD_BUILD_DIR} && \
    cmake ${SCD_SRC_DIR}
RUN ln -s -f /usr/lib64/libGL.so.1 /usr/lib64/libGL.so && \
    ln -s -f /usr/lib64/libexpat.so.1 /usr/lib64/libexpat.so && \
    for libTag in Xmu ICE Xext X11 Xt SM; do ln -s -f /usr/lib64/lib${libTag}.so.6 /usr/lib64/lib${libTag}.so; done
RUN cd ${SCD_BUILD_DIR} && \
    make -j4
