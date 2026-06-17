# JSON Likelihoods for ATLAS SUSY displaced leptons search (SUSY-2018-14)

The JSON likelihoods are serialized for each signal region: [SRee](SRee), [SRmm](SRmm), [Stau](Stau), as well as for the superpartners of the left/right handed leptons, denoted as, for example, `[region]_left`, and `[region]_right`. This is done by providing a background-only workspace containing the signal channels for each region at `[region]_BkgOnly.json` as well as patch files for each mass point on the signal phase-space explored in the analysis, at `[region]_patchset.json`.

Each [jsonpatch](http://jsonpatch.com/) file follows the format `DisplacedLeptons_[region]_[signal]_[m]_[tau].json` where `signal` is either `Slepton` (for selectron and smuon) or `Stau`, `m` is the mass of the slepton, and `tau` is the lifetime in picoseconds.

One can list all the patches available in a patchset.json file using for example [jq](https://stedolan.github.io/jq/), like:
```
jq '.patches[].metadata.name' SRee_patchset.json
```
shown here for SRee only.

## Producing signal workspaces

Signal workspaces can be produced using [python jsonpatch](https://python-json-patch.readthedocs.io/en/latest/) and [pyhf](https://scikit-hep.org/pyhf/) like:

```
jsonpatch SRee_bkgonly.json <(pyhf patchset extract SRee_patchset.json --name "DisplacedLeptons_SRee_Slep_1000_10") > SRee_DisplacedLeptons_SRee_Slep_1000_10.json
```
here shown as an example for a signal model `DisplacedLeptons_SRee_Slep_1000_10`.

## Computing signal workspaces

For example, with [pyhf](https://scikit-hep.org/pyhf/), one can do any of the following:

```
pyhf cls SRee_bkgonly.json -p <(pyhf patchset extract SRee_patchset.json --name "DisplacedLeptons_SRee_Slep_1000_10") 

jsonpatch SRee_bkgonly.json <(pyhf patchset extract SRee_patchset.json --name "DisplacedLeptons_SRee_Slep_1000_10") | pyhf cls

pyhf cls SRee_DisplacedLeptons_SRee_Slep_1000_10.json
```
