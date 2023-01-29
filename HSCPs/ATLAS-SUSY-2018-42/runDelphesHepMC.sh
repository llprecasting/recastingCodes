#!/bin/sh

currentDIR="$( pwd )"
delphesDIR=$currentDIR/DelphesHSCP
pythiaDIR=$currentDIR/MG5/HEPTools/pythia8


Help()
{
   # Display Help
   echo
   echo "Syntax: ./runDelphesHepMC <delphes card> <output file> <input file>"
   echo
}


#Make sure pythia can be found by Delphes
export LD_LIBRARY_PATH=$pythiaDIR/lib:$LD_LIBRARY_PATH
#Make sure Delphes can be found by ROOT
export ROOT_INCLUDE_PATH=$delphesDIR/external

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then 
  echo "Export the necessary paths to ROOT_INCLUDE_PATH and LD_LIBRARY_PATH and run $delphesDIR/DelphesHepMC2."
  Help;
  exit
fi


if [ $# -ne "3" ]; then 
  Help;
  exit
fi



if [ ! -f "$currentDIR/$1" ]; then
  echo "Delphes card $currentDIR/$1 DOES NOT exists.\n Please, pass a valid path as the first argument";
  exit
fi

if [ -f "$currentDIR/$2" ]; then 
#  echo "Output file $currentDIR/$2 already exists. DELETE file?[y/n]";
#  read answer;
#  if echo "$answer" | grep -iq "^y" ;then
    rm $currentDIR/$2;
#  else
#    exit;    
#  fi
fi

if [ ! -f "$currentDIR/$3" ]; then 
  echo "Input file $currentDIR/$3 DOES NOT exists.\n Please, pass a valid path as the third argument";
  exit;    
fi




cd $delphesDIR
./DelphesHepMC2 $currentDIR/$1 $currentDIR/$2 $currentDIR/$3
cd $currentDIR
echo "Output written to $currentDIR/$2"
