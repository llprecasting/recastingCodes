To process the recasting, please follow these instructions.

- First of all, you will need to install this version of madgraph : MG5_aMC_v3.4.2.tar.gz from : https://launchpad.net/mg5amcnlo/+download

- The required UFO files are available in the /UFO folder 

- install `lhapdf` and `pythia8` as:
```
./bin/mg5_aMC
MG5_aMC> install lhapdf6
MG5_aMC> install pythia8
```

- then you can use the scripts to generate and process events like so, depending on the model you wish to choose.

```
python3 Writing_Scripts_MG+P8_Single_HS1S2.py gg 400_100_100_100.0_100.0 20000 <path to madgraph> <path to where to store the HEPMC files> 100.0_100.0
cd /users/divers/atlas/corpe/scratch/RecastingLouisePaper
/AtlasDisk/user/corpe/LLPRecasting2024/MG5_aMC_v3_4_2/bin/mg5_aMC -f /users/divers/atlas/corpe/scratch/RecastingLouisePaper/script_ggH400_S1100_S2100_gg_100.0_100.0.txt
cd /AtlasDisk/user/corpe/LLPRecasting2024/MG5_aMC_v3_4_2/recastingCodes/CalRatioDisplacedJet/
python3 /AtlasDisk/user/corpe/LLPRecasting2024/MG5_aMC_v3_4_2/recastingCodes/CalRatioDisplacedJet//Computation_Map_Single_ALL.py gg 400_100_100_100.0_100.0 HS1S2 /AtlasDisk/user/corpe/LLPRecasting2024/MG5_aMC_v3_4_2 /users/divers/atlas/corpe/scratch/RecastingLouisePaper 100.0_100.0

python3 Writing_Scripts_MG+P8_Single_HSS.py gg 400_100 20000 <path to madgraph> <path to where to store the HEPMC files>
cd <path to where to store the HEPMC files>
<path to madgraph>/bin/mg5_aMC -f <path to where to store the HEPMC files>/script_ggH400_S100_gg.txt
cd <here>
python3 Computation_Map_Single_ALL.py gg 400_100 HSS <path to madgraph> <path to where to store the HEPMC files>
```
