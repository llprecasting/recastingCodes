#!/bin/sh

homeDIR="$( pwd )"


#Get HepMC tarball
hepmc="hepmc2.06.11.tgz"
echo -n "Install HepMC2 (y/n)? "
read answer
if echo "$answer" | grep -iq "^y" ;then
	mkdir hepMC_tmp
	URL=http://hepmc.web.cern.ch/hepmc/releases/$hepmc
	echo "[installer] getting HepMC"; wget $URL 2>/dev/null || curl -O $URL; tar -zxf $hepmc -C hepMC_tmp;
	mkdir HepMC-2.06.11; mkdir HepMC-2.06.11/build; mkdir HepMC2;
	echo "Installing HepMC in ./HepMC";
	cd HepMC-2.06.11/build;
	../../hepMC_tmp/HepMC-2.06.11/configure --prefix=$homeDIR/HepMC2 --with-momentum=GEV --with-length=MM;
	make;
	make check;
	make install;

	#Clean up
	cd $homeDIR;
	rm -rf hepMC_tmp; rm $hepmc; rm -rf HepMC-2.06.11;
fi


#Get FastJet tarball
fastjet="fastjet-3.4.2.tar.gz"
URL=http://fastjet.fr/repo/$fastjet
echo -n "Install FastJet (y/n)? "
read answer
if echo "$answer" | grep -iq "^y" ;then
	mkdir fastjet_tmp
	echo "[installer] getting Fastjet"; wget $URL 2>/dev/null || curl -O $URL; tar -zxf $fastjet -C fastjet_tmp;
	mkdir fastjet-3.4.2;
	echo "Installing Fastjet in fastjet-3.2.2";
	cd fastjet_tmp/fastjet-3.4.2;
	./configure  --enable-allplugins --prefix=$homeDIR/fastjet-3.4.2/;
	make install;

	#Clean up
	cd $homeDIR
	rm -rf fastjet_tmp; rm $fastjet;
fi

#Get pythia tarball
pythia="pythia8310.tgz"
URL=https://pythia.org/download/pythia83/$pythia
echo -n "Install Pythia (y/n)? "
read answer
if echo "$answer" | grep -iq "^y" ;then
	if hash gzip 2>/dev/null; then
		mkdir pythia8;
		echo "[installer] getting Pythia"; wget $URL 2>/dev/null || curl -O $URL; 
		tar -zxf $pythia -C pythia8 --strip-components 1;
		echo "Installing Pythia in pythia8";
		cd pythia8;
		./configure --with-root=$ROOTSYS --prefix=$homeDIR/pythia8 --with-gzip --with-fastjet3=../fastjet-3.4.2 --with-hepmc2=$homeDIR/HepMC2
		make -j4; make install;
		cd $homeDIR
		rm $pythia;
	else
		echo "[installer] gzip is required. Try to install it with sudo apt-get install gzip";
	fi
fi



