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
RUN  dnf -y copr enable averbyts/HEPrpms && yum -y  install geant4-devel geant4 fastjet fastjet-devel \
     	    	 			    	    clhep clhep-devel PTL PTL-devel \
						    pythia8-devel pythia8 python3-lhapdf lhapdf lhapdf-devel \
						    HepMC3* HepMC HepMC-devel \
						    &&  yum -y clean all

ENV Geant4_DataDir=/usr/share/Geant4/data
RUN mkdir ${Geant4_DataDir}

RUN cd ${Geant4_DataDir} && \
    wget http://cern.ch/geant4-data/datasets/G4NDL.4.7.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4EMLOW.8.5.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4PII.1.3.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4PhotonEvaporation.5.7.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4RadioactiveDecay.5.6.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4RealSurface.2.2.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4SAIDDATA.2.0.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4ABLA.3.3.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4INCL.1.2.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4ENSDFSTATE.2.3.tar.gz && \
    wget http://cern.ch/geant4-data/datasets/G4PARTICLEXS.4.0.tar.gz && \
    tar -xzf G4NDL.4.7.tar.gz && \
    tar -xzf G4EMLOW.8.5.tar.gz && \
    tar -xzf G4PII.1.3.tar.gz && \
    tar -xzf G4PhotonEvaporation.5.7.tar.gz && \
    tar -xzf G4RadioactiveDecay.5.6.tar.gz && \
    tar -xzf G4RealSurface.2.2.tar.gz && \
    tar -xzf G4SAIDDATA.2.0.tar.gz && \
    tar -xzf G4ABLA.3.3.tar.gz && \
    tar -xzf G4INCL.1.2.tar.gz && \
    tar -xzf G4ENSDFSTATE.2.3.tar.gz && \
    tar -xzf G4PARTICLEXS.4.0.tar.gz && \
    rm -f *tar.gz

ENV G4ABLADATA=${Geant4_DataDir}/G4ABLA3.3
ENV G4INCLDATA=${Geant4_DataDir}/G4INCL1.2
ENV G4LEDATA=${Geant4_DataDir}/G4EMLOW8.5
ENV G4PARTICLEXSDATA=${Geant4_DataDir}/G4PARTICLEXS4.0
ENV G4REALSURFACEDATA=${Geant4_DataDir}/RealSurface2.2
ENV G4LEVELGAMMADATA=${Geant4_DataDir}/G4PhotonEvaporation5.7
ENV G4PIIDATA=${Geant4_DataDir}/G4PII1.3
ENV G4SAIDXSDATA=${Geant4_DataDir}/G4SAIDDATA2.0
ENV G4ENSDFSTATEDATA=${Geant4_DataDir}/G4ENSDFSTATE2.3
ENV G4NEUTRONHPDATA=${Geant4_DataDir}/G4NDL4.7
ENV G4RADIOACTIVEDATA=${Geant4_DataDir}/G4RadioactiveDecay5.6

ENV WKDIR=/work
RUN mkdir ${WKDIR}
WORKDIR ${WKDIR}

COPY . ${WKDIR}
ENV COCOA_BASE_DIR=${WKDIR}/COCOA
ENV COCOA_BUILD_DIR=${WKDIR}/cocoa-build
RUN mkdir ${COCOA_BUILD_DIR}
RUN cd ${COCOA_BUILD_DIR} && \
    cmake ${COCOA_BASE_DIR}
RUN cd ${COCOA_BUILD_DIR} && \
    make -j10
