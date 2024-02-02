# Validation #

## Running ##

For obtaining the validation results shown below, the following steps should be taken:

 1. Generate MC event samples for the strong or electroweak benchmarks. This can be done running the the script [runScanMG5.py](../runScanMG5.py) with the [strong scan parameters](./scan_parameters_strong.ini) or the [ewk scan parameters](./scan_parameters_ewk.ini).
 
 2. Compute the cut flows for each benchmark running:

```
atlas_susy_2018_13_CutFlow.py -f <Delphes root_file>  -m <model> -s <signal region>
```

where model should be "strong" or "ewk" and signal region = "HighPT" or "Trackless". The cutflow output will be printed to the screen and stored in a pickle file in the form of a Pandas DataFrame.

3. The cut flows for each benchmark can be combined running:

```
atlas_susy_2018_13_CombineData.py -f <list of pickle files> -o <output file>
```

4. Finally, the Jupyter notebook [getCutFlowTables.ipynb](./getCutFlowTables.ipynb) can be used to read the cut flows and compare to the official ATLAS values.

## Validation Results ##


### Cut-flow Comparison


 * Strong production: values correspond to Recast (ATLAS) acceptance and final acceptance times efficiency.


|                              |               |               |               |               |
|:-----------------------------|:--------------|:--------------|:--------------|:--------------|
| $m_{\tilde g} (GeV)$         | 2000.0        | 2000.0        | 2400.0        | 2000.0        |
| $m_{\tilde \chi_1^0} (GeV)$  | 850.0         | 50.0          | 200.0         | 1250.0        |
| $\tau(\tilde \chi_1^0) (ns)$ | 0.01          | 0.1           | 1.0           | 10.0          |
| Total                        | 1.0           | 1.0           | 1.0           | 1.0           |
| Jet selection                | 0.999 (0.999) | 0.967 (0.966) | 0.985 (0.972) | 0.999 (0.961) |
| $R_{xy},z <$ 300 mm          | 0.999 (0.999) | 0.788 (0.787) | 0.442 (0.447) | 0.311 (0.317) |
| $R_{DV} > 4$ mm              | 0.298 (0.296) | 0.772 (0.770) | 0.434 (0.438) | 0.304 (0.309) |
| $d_0 > 2$ mm                 | 0.296 (0.296) | 0.767 (0.756) | 0.433 (0.437) | 0.304 (0.309) |
| $nTracks >= 5$               | 0.296 (0.296) | 0.766 (0.755) | 0.433 (0.437) | 0.304 (0.309) |
| mDV > 10 GeV                 | 0.296 (0.296) | 0.760 (0.747) | 0.433 (0.437) | 0.304 (0.309) |
| final Acc*Eff                | 0.274 (0.278) | 0.140 (0.144) | 0.112 (0.115) | 0.090 (0.092) |


 * Electroweak production: values correspond to Recast (ATLAS) acceptance and final acceptance times efficiency.

|                              |               |               |               |               |
|:-----------------------------|:--------------|:--------------|:--------------|:--------------|
| $m_{\tilde \chi_1^0} (GeV)$  | 500.0         | 500.0         | 1300.0        | 1300.0        |
| $\tau(\tilde \chi_1^0) (ns)$ | 0.1           | 1.0           | 0.1           | 1.0           |
| Total                        | 1.0           | 1.0           | 1.0           | 1.0           |
| Jet selection                | 0.502 (0.495) | 0.513 (0.501) | 0.984 (0.968) | 0.986 (0.985) |
| $R_{xy},z <$ 300 mm          | 0.502 (0.495) | 0.416 (0.410) | 0.984 (0.968) | 0.924 (0.921) |
| $R_{DV} > 4$ mm              | 0.465 (0.465) | 0.403 (0.398) | 0.865 (0.859) | 0.905 (0.899) |
| $d_0 > 2$ mm                 | 0.464 (0.465) | 0.402 (0.398) | 0.863 (0.859) | 0.905 (0.899) |
| $nTracks >= 5$               | 0.464 (0.465) | 0.402 (0.398) | 0.863 (0.859) | 0.905 (0.899) |
| mDV > 10 GeV                 | 0.464 (0.465) | 0.402 (0.398) | 0.863 (0.859) | 0.905 (0.899) |
| final Acc*Eff                | 0.288 (0.311) | 0.142 (0.143) | 0.104 (0.122) | 0.080 (0.083) |