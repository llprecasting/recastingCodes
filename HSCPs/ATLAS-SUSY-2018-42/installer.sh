#!/bin/sh

homeDIR="$( pwd )"

echo "Installation will take place in $homeDIR"

echo "[Checking system dependencies]"
PKG_OK=$(dpkg-query -W -f='${Status}' autoconf 2>/dev/null | grep -c "ok installed")
if test $PKG_OK = "0" ; then
  echo "autoconf not found. Install it with sudo apt-get install autoconf."
  exit
fi
PKG_OK=$(dpkg-query -W -f='${Status}' libtool 2>/dev/null | grep -c "ok installed")
if test $PKG_OK = "0" ; then
  echo "libtool not found. Install it with sudo apt-get install libtool."
  exit
fi
PKG_OK=$(dpkg-query -W -f='${Status}' gzip 2>/dev/null | grep -c "ok installed")
if test $PKG_OK = "0" ; then
  echo "gzip not found. Install it with sudo apt-get install gzip."
  exit
fi
PKG_OK=$(dpkg-query -W -f='${Status}' bzr 2>/dev/null | grep -c "ok installed")
if test $PKG_OK = "0" ; then
  echo "bzr not found. Install it with sudo apt-get install bzr."
  exit
fi

cd $homeDIR

madgraph="MG5_aMC_v3.4.2.tar.gz"
URL=https://launchpad.net/mg5amcnlo/3.0/3.4.x/+download/$madgraph
echo -n "Install MadGraph (y/n)? "
read answer
if echo "$answer" | grep -iq "^y" ;then
	mkdir MG5;
	echo "[installer] getting MadGraph5"; wget $URL 2>/dev/null || curl -O $URL; tar -zxf $madgraph -C MG5 --strip-components 1;
	cd ./MG5/bin;
	echo "[installer] installing HepMC under MadGraph5"
	echo "install hepmc\nexit\n" > mad_install.txt
	./mg5_aMC -f mad_install.txt
	cd ../;
	"[installer] Trying to install lhapdf 6.5.3. under MadGrapg5";
	sed -i "s/'version':       '6.3.0'/'version':       '6.5.3'/g" HEPTools/HEPToolsInstallers/HEPToolInstaller.py;
	python3 ./HEPTools/HEPToolsInstallers/HEPToolInstaller.py lhapdf6;
	if [ ! -f "./HEPTools/lhapdf6_py3/bin/lhapdf-config" ]; then	
	    echo "LHAPDF6 installation failed. Try to install it manually";
	    exit;
	fi
	cd bin;
        echo "[installer] installing Pythia8 under MadGraph5";
	echo "install pythia8\nexit\n" > mad_install.txt;
	./mg5_aMC -f mad_install.txt;
	rm mad_install.txt;
	cd $homeDIR;
	sed  "s|homeDIR|$homeDIR|g" mg5_configuration.txt > ./MG5/input/mg5_configuration.txt;
        rm $madgraph;
fi

#Get HepMC tarball
#hepmc="hepmc2.06.11.tgz"
#echo -n "Install HepMC2 (y/n)? "
#read answer
#if echo "$answer" | grep -iq "^y" ;then
#	mkdir hepMC_tmp
#	URL=http://hepmc.web.cern.ch/hepmc/releases/$hepmc
#	echo "[installer] getting HepMC"; wget $URL 2>/dev/null || curl -O $URL; tar -zxf $hepmc -C hepMC_tmp;
#	mkdir HepMC-2.06.11; mkdir HepMC-2.06.11/build; mkdir HepMC2;
#	echo "Installing HepMC in ./HepMC";
#	cd HepMC-2.06.11/build;
#	../../hepMC_tmp/HepMC-2.06.11/configure --prefix=$homeDIR/HepMC2 --with-momentum=GEV --with-length=MM;
#	make;
#	make check;
#	make install;
#
#	#Clean up
#	cd $homeDIR;
#	rm -rf hepMC_tmp; rm $hepmc; rm -rf HepMC-2.06.11;
#fi

#Get pythia tarball
#pythia="pythia8308.tgz"
#URL=https://pythia.org/download/pythia83/$pythia
#echo -n "Install Pythia (y/n)? "
#read answer
#if echo "$answer" | grep -iq "^y" ;then
#	if hash gzip 2>/dev/null; then
#		mkdir pythia8;
#		echo "[installer] getting Pythia"; wget $URL 2>/dev/null || curl -O $URL; tar -zxf $pythia -C pythia8 --strip-components 1;
#		echo "Installing Pythia in pythia8";
#		cd pythia8;
#		./configure --with-hepmc2=$homeDIR/HepMC2 --with-root=$ROOTSYS --prefix=$homeDIR/pythia8 --with-gzip
#		make -j4; make install;
#		cd $homeDIR
#		rm $pythia;
#	else
#		echo "[installer] gzip is required. Try to install it with sudo apt-get install gzip";
#	fi
#fi


echo -n "Install Delphes (y/n)? "
repo=https://github.com/delphes/delphes
URL=http://cp3.irmp.ucl.ac.be/downloads/$delphes
read answer
if echo "$answer" | grep -iq "^y" ;then
#  latest=`git ls-remote --sort="version:refname" --tags $repo  | grep -v -e "pre" | grep -v -e "\{\}" | cut -d/ -f3- | tail -n1`
#  echo "[installer] Cloning Delphes version $latest";
#  git clone --branch $latest https://github.com/delphes/delphes.git Delphes
#  echo "[installer] Adding HSCP module to Delphes";
#  cp delphesHSCP/Makefile Delphes/;
#  cp delphesHSCP/HSCPFilter.cc delphesHSCP/HSCPFilter.h delphesHSCP/ModulesLinkDef.h Delphes/modules/;
#  cd Delphes;
#  echo "[installer] installing Delphes";

# Check if pythia8 has been installed under MadGraph5
  pythiaDir=$homeDIR/MG5/HEPTools/pythia8 
  if [ ! -d "$pythiaDir" ]; then
    echo "Delphes should be installed after hepmc, lhapdf6 and pythia8 were installed in MadGraph."
    exit
  fi
  echo "[installer] Installing DelphesHSCP";    
  tar -zxf DelphesHSCP.tar.gz;
  cd DelphesHSCP;
  export PYTHIA8=$pythiaDir;
  make HAS_PYTHIA8=true;
  rm -rf .git
  cd $homeDIR;
fi

echo "\n[installer] For running Delphes the following env variables should be set:\n\n export LD_LIBRARY_PATH=$homeDIR/MG5/HEPTools/pythia8/lib"
echo "\nand for reading Delphes produced ROOT files:\n\n export ROOT_INCLUDE_PATH=$homeDIR/DelphesHSCP/external\n"
echo "\n\n or run source setenv.sh\n\n"

cd $homeDIR
