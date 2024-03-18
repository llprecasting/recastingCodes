# Validation #

## Running ##

For obtaining the validation results shown below, the following steps should be taken:

 1. Generate MC event samples for the strong or electroweak benchmarks. This can be done running the the script [runScanMG5.py](../runScanMG5.py) with the [gluino scan parameters](./scan_parameters_gluino.ini).
 
 2. Compute the signal efficiencies and yield:

```
atlas_susy_2018_13_Recast.py -f <Delphes root_file>
```

The recasting output will be printed to the screen and stored in a pickle file in the form of a Pandas DataFrame.

3. The efficiencies and signal yields for each benchmark can be combined running:

```
atlas_susy_2016_13_CombineData.py -f <list of pickle files> -o <output file>
```

4. Finally, the Jupyter notebooks [validation-mlsp100.ipynb](./validation-mlsp100.ipynb) and  [validation-mlsp1300.ipynb](./validation-mlsp1300.ipynb) can be used to produce the
results below

## Validation Results ##


The validation of the signal efficiencies (efficiency times acceptance)
for direct production of staus can be found in the [validation folder](validation).
The output (.eff files) was generated using the SLHA files and pythia8.cfg and parameters.ini files stored in the folder.
For instance, the following validation plot can be generated running this [ipython notebook](validation/validation-mlsp100.ipynb):

 * mgluino = 1.4 TeV, mLSP = 100 GeV:

![Alt text](validationPlot_mlsp100.png?raw=true "Validation Plot")


* mgluino = 1.4 TeV, mLSP = 1300 GeV:

![Alt text](validationPlot_dm100.png?raw=true "Validation Plot")



### Cut-flow Comparison


|             | ATLAS  | Recasting  |   ATLAS | Recasting  |
| ----------- | ------ | ---------- | ------- | ---------- |
|  **mGluino, mLSP, tau**             | **2000,100,1ns**  | **2000,100,1ns**  |   **2000,1800,1ns** | **2000,1800,1ns**  |
| Initial Events                 |   32   |   32   |   32   |   32  |
| Trigger-based data reduction   |   32   |        |   27   |       |
| Event cleaning                 |   32   |        |   27   |       |
| Good Runs List   |   32   |   --   |   27   |   --  |
| Primary vertex   |   32   |   --   |   27   |   --  |
| NCB veto   |   32   |   --   |   26   |   --  |
| MET trigger   |   31   |   --   |   24   |   --  |
| MET filter (MET > 200 GeV for recast)   |   31   |   30.4   |   17   |   15  |
| Offline MET (event selection eff. for recast)   |   29   |   26.8   |   7   |  9.3  |
|  |
| DV fiducial acceptance   |   28   |   --   |   6   |   --  |
| DV fit quality   |   27   |   --   |   6   |   --  |
| DV displacement   |   27   |   --   |   6   |   --  |
| Material veto   |   22   |   --   |   5   |   --  |
| Disabled module veto   |   22   |   --   |   5   |   --  |
| DV track multplicity   |   15   |   --   |   3   |   --  |
| DV mass (vertex level eff. for recast)   |   14   |   14.4   |   2   |   3.6  | 	
