FROM gitlab-registry.cern.ch/linuxsupport/rpmci/builder-cs9

RUN  yum -y  install  dnf-plugins-core --exclude=*uploa* --exclude=*product* --exclude=*subscr* epel*
RUN  yum -y  install  bc make cmake binutils git wget diffutils file sed gawk grep which xerces-c xerces-c-devel \
     	     	      emacs mesa-libGL-devel libXmu-devel expat expat-devel \
                      gcc-gfortran gcc-c++ bzip2 openssl-devel openssl jsoncpp jsoncpp-devel \
                      libzip-devel  zlib zlib-devel \
                      root root-* \
                      python3 python3-devel &&  yum -y  clean all
RUN ln -s /usr/include/boost/json /usr/include/json

#Packages coming from HEPrpms ->
RUN  dnf -y copr enable averbyts/HEPrpms && yum -y  install fastjet fastjet-devel \
     	    	 			    	    clhep clhep-devel PTL PTL-devel \
						    pythia8-devel pythia8 python3-lhapdf lhapdf lhapdf-devel \
						    HepMC3* HepMC HepMC-devel \
						    &&  yum -y clean all



ENV WKDIR=/work
RUN mkdir ${WKDIR}
WORKDIR ${WKDIR}

###########
# Geant 4 #
###########

RUN cd ${WKDIR}
RUN wget http://cern.ch/geant4-data/releases/lib4.10.7.p03/Linux-g++8.3.0-CC7.tar.gz
RUN tar -xzf Linux-g++8.3.0-CC7.tar.gz

ENV Geant4_DIR=${WKDIR}/Geant4-10.7.3-Linux
ENV Geant4_InstallBaseDir=${Geant4_DIR}/share/Geant4-10.7.3
ENV Geant4_InstallDir=${Geant4_DIR}/geant4-10.7.3
RUN mkdir ${Geant4_InstallDir}
RUN wget http://cern.ch/geant4-data/releases/geant4.10.07.p03.tar.gz
RUN tar -xzf geant4.10.07.p03.tar.gz

ENV G4_BuildDir=${WKDIR}/g4_build
RUN mkdir ${G4_BuildDir} && \
    cd ${G4_BuildDir} && \
    cmake -DCMAKE_INSTALL_PREFIX:PATH=${Geant4_InstallDir=$} -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_INSTALL_DATA=ON ${WKDIR}/geant4.10.07.p03 && \
    make -j16 && \
    make install

ENV Geant4_DataDir=/work/Geant4-10.7.3-Linux/geant4-10.7.3/share/Geant4-10.7.3/data
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

#######
# COCOA #
#######

ENV HEPMC_HOME=/usr

ENV COCOA_BASE_DIR=${WKDIR}/COCOA
ENV COCOA_SRC_DIR=${COCOA_BASE_DIR}/COCOA
ENV COCOA_BUILD_DIR=${COCOA_BASE_DIR}/build
RUN mkdir ${COCOA_BASE_DIR} ${COCOA_BUILD_DIR}

COPY . ${COCOA_BASE_DIR}
RUN cd ${COCOA_BUILD_DIR} && \
    cmake ${COCOA_SRC_DIR}
RUN cd ${COCOA_BUILD_DIR} && \
    make -j16
