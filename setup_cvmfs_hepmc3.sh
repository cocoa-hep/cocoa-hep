export CURRENTDIR=$(pwd)
source /cvmfs/sft.cern.ch/lcg/releases/LCG_96b/Geant4/10.05.p01.1/x86_64-centos7-gcc8-opt/Geant4-env.sh
cd /cvmfs/sft.cern.ch/lcg/releases/LCG_96b/Geant4/10.05.p01.1/x86_64-centos7-gcc8-opt/share/Geant4-10.5.1/geant4make/
source /cvmfs/sft.cern.ch/lcg/releases/LCG_96b/Geant4/10.05.p01.1/x86_64-centos7-gcc8-opt/share/Geant4-10.5.1/geant4make/geant4make.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_96b/ROOT/6.18.04/x86_64-centos7-gcc8-opt/ROOT-env.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_96b/jsoncpp/1.8.4/x86_64-centos7-gcc8-opt/jsoncpp-env.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_96b/CMake/3.14.3/x86_64-centos7-gcc8-opt/CMake-env.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_96b/cmaketools/1.8/x86_64-centos7-gcc8-opt/cmaketools-env.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_96b/fastjet/3.3.2/x86_64-centos7-gcc8-opt/fastjet-env.sh
#source /cvmfs/sft.cern.ch/lcg/releases/LCG_96b/MCGenerators/pythia8/301/x86_64-centos7-gcc8-opt/pythia8env-genser.sh
# < install hepmc3 if not done yet >
export HEPMC_HOME=/home/pr2329/cocoa-hep/hepmc3_noPythia8Interface
ln -s ${HEPMC_HOME}/lib64 ${HEPMC_HOME}/lib
# < install pythia 8, using './configure --with-hepmc3=${HEPMC_HOME}' >
export PYTHIA8_HOME=/home/pr2329/pythia8311
export PYTHIA8DATA=$PYTHIA8_HOME/share/Pythia8/xmldoc   
cd $CURRENTDIR
