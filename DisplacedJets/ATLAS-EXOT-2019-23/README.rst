To process the recasting, please follow these instructions.

- First of all, you will need to install this version of madgraph : MG5_aMC_v3.4.2.tar.gz from : https://launchpad.net/mg5amcnlo/+download

- In the same folder, download and place the model named : HAHM_MG5model_v3 from http://insti.physics.sunysb.edu/~curtin/hahm_mg.html. Extract the tarball in the MG5_aMC_v3_4_2 folder.

- You will also need the folder 'ATLAS_data' in this github repository and place this folder in your MG5_aMC_v3_4_2 folder. It contains files from HEP_Data from the ATLAS public results.

- Launch madgraph by writing, in a konsole (launched in the MG5_aMC_v3_4_2 folder) : ./bin/mg5_aMC

- After entering in the MG konsole, install LHAPDF6 and Pythia8 by writing :

MG5_aMC> install lhapdf6

And after :

MG5_aMC> install pythia8

It will take some minutes.

- You will also need to convert the model file by writing :

MG5_aMC> convert model ./HAHM_MG5model_v3/HAHM_variableMW_v3_UFO

MG is ready to generate events on your computer ! To continu, please be sure that everything that you download from the GitHub repository should be IN the folder MG5_aMC_v3_4_2 ( where is your madgraph )

- Now, you can generate the scripts that MG will read by launching the function Writing_Scripts_MG+P8.py in the konsole with :

python3 Writing_Scripts_MG+P8.py

It will write five scripts and launch the associated madgraph generation. This generation is long. After approximately 1h40-2h of running, all files should be ready.

- To process the data and obtain plots and limits, you can now launch in the konsole :

python3 Computation_Map.py

It will also be long ( more than 2h ) and at the end, you will have the plots, the limits and the values in a text files, so that you do not have to redo the entire run.

