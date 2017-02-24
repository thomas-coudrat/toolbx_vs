## Synopsis

**toolbx_vs** is a set of tools that help setup, run, monitor and analyze
virtual screens (VS) on a cluster.

## Examples

Below are examples that showcase each of the functionality of toolbx_vs. Each
script presented here has a --help flag that can be used to get more
information on its use.

### Preparation

**Create chemical library index file**
Create a .ndx file with suffix "\_clusterA" on the small molecule library
chemical_lib.sdf
```
vs_index.py chemical_lib.sdf clusterA
```

**Create maps of binding pocket for docking**
Create maps for docking of the protein receptor.ob to be screened by the
chemical library chemical_lib_clusterA.inx. The database type is set to 3D
structure, and the maps are created using "ligand mode" which generates the maps
around the bound ligand in receptor.ob (and remove that ligand).
```
vs_maps.py receptor.ob chemical_lib_clusterA.inx 3D ligand
```

**Setup virtual screen parameters**
Builds the VS to screen molecules 200 to 1000 of the chemical library, splitting
it into slices of 100 (8 slices). The number of repeats is set to 3 and
thoroughness to 10. The time limit on the cluster is set to 24 hours, and the
SLURM scheduling system is to be used. Finally vs_setup/ is the name of the
directory that contains the receptor, its maps and the chemical library .inx
files.
```
vs_build.py 200 1000 100 3 10. 0-24:00:00 vs_setup slurm
```

### Execution

**Execute virtual screen on a cluster**
Run the VS that was setup above by submitting the name of the directory where
it was built, and specifying the scheduling system (SLURM in this case).
```
vs_submit.py my_vs_experiment/ slurm
```

**Print report on virtual screen progress**
Print a report of the process of the VS on the cluster. Run in a VS directory.
```
vs_report.py
```

### Analysis

**Extract virtual scree results**
Extract and consolidate the VS ranked results from the multiple parallel docking
runs that were completed.
```
vs_results.py my_vs_experiment/
```

**Plot ROC curve**
This plots a ROC curve molecules 200 to 600 as true positives and 601 to 1000 as
false positives. The figure is in log scale, the curve is red and continuous.
```
vs_plot_roc.py 'My VS experiment' 'Receptor' 'my_vs_experiment/results_receptor.csv' 'red' 'cont' '200-600' '601-1000' '% receptor agonists' '% receptor inhibitors' -log
```

**Plot enrichment factor bargraph**
Plot a enrichment factor for defined sets of ligands in "ligand_types.json".
```
vs_plot_ef.py 'My VS experiment' 'Receptor' 'my_vs_experiment/results_receptor.csv' 'red' '200-600' '601-1000' ligand_types.json
```

**Extract docked poses**
Extract docked poses from the VS results for analysis. This extracts the top 10
poses and poses for ligands with IDs 315 and 2017.
```
vs_poses.py 'my_vs_experiment/results_receptor.csv' 10 --ligIDs 315,2017
```

## Motivation
This set of tools simplifies the management of a VS on an HPC cluster with ICM.
These were created for VS experiments performed by Thomas Coudrat during his
PhD. Features were added as needed for the study. The scripts don't provide
access to all VS functionalities available from the ICM docking software.

## Contributions
Feel free to create an issue or submit a pull request. You can also contact me
(Thomas Coudrat) if you have questions related to this project.

## Installation
* Install Anaconda for Python 3.5
* Install ICM 3.8-4
    * set environment variable: export ICMHOME="path/to/icmDirectory"
    * Requires software license

## License
This project is licensed under the MIT license
