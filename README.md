# DFMFON
[![Last Commit](https://img.shields.io/github/last-commit/smbeselly/DFMFON/commits/main)](
https://github.com/smbeselly/DFMFON/commits/main)
[![GitHub release](https://img.shields.io/github/release/smbeselly/DFMFON)](https://GitHub.com/smbeselly/DFMFON/releases/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

DFMFON, stands for Delft3D-Flexible Mesh (DFM) and MesoFON (MFON) is an open-source software written in Python to mechanistically simulate the Mangrove and Hydromorphology development. We achieve that by coupling multi-paradigm of individual-based mangrove model [MFON](http://mesofon.org/index.php) and process-based hydromorphodynamic model [DFM](https://oss.deltares.nl/web/delft3dfm).

This model is developed as part of the PhD research by Sebrian Beselly at [IHE Delft](https://www.un-ihe.org/department/coastal-and-urban-risk-resilience) and [Delft University of Technology](https://www.tudelft.nl/).

Publication describing the DFMFON

[Modelling mangrove-mudflat dynamics with a coupled individual-based-hydro-morphodynamic model](https://doi.org/10.1016/j.envsoft.2023.105814) (open access)

# What is DFMFON

As climate-change-driven extremes potentially make coastal areas more vulnerable, mangroves can help sustainably protect the coasts. There is a substantial understanding of both mangrove dynamics and hydro-morphodynamic processes. However, the knowledge of complex eco-geomorphic interactions with physical-environmental stressors remains lacking. We introduce a novel coupled modelling approach consisting of an individual-based mangrove (mesoFON) and a process-based hydromorphodynamic model (Delft3D-FM). This coupled model is unique because it resolves spatiotemporal processes, including tidal, seasonal, and decadal environmental changes (water level, flow, sediment availability, and salinity) with full life-stages (propagule, seedling, sapling, mature) mangrove interaction. It allows us to mechanistically simulate forest expansion, retreat, and colonisation influenced by and with feedback on physical-environmental drivers. The model is applied in a schematized mixed fluvial-tidal deltaic mangrove forest in dominantly muddy sediment inspired by the prograding delta of Porong, Indonesia. Model results successfully reproduce observed mangrove extent development, age-height relationship, and morphodynamic delta features.

# 1. Installation
## 1.1. Install Git
It is wise to have the git installed in your PC.

You can directly download the executable in WIindows system through this [link](https://git-scm.com/download/win), or refer to the [Git official website](https://git-scm.com/) for other information. 
## 1.2. Create a virtual environment with Anaconda
To run the model, it is advised to create a virtual environment with Anaconda. We have prepared the explicit text file containing the spec list.

First, clone the main repository or download the ZIP file.

Change the working directory to the repo location.

In the current working directory, run:

    conda create --name dfmfon --file spec-virtual_env.txt

Once it is finished, activate the newly created virtual environment and install the BMI Python for Delft3D FLexible Mesh

    conda activate dfmfon
    pip install bmi-python

Now all the required packages have been installed and you are all set to run the DFMFON model. Always make sure to activate the `DFMFON` environment.

# 2. Usage
An example to run DFMFON is provided in the `example' folder.

We provide *run_example* to try the tool and **folder** template, which you can use for your repeatable runs.

## 2.1. Preparation Run
DFMFON currently has two main script, `01_prep_run.py` and `02_DFMFON_run.py`.

At first, open `01_prep_run.py` in your favorite IDE, I suggest [Spyder](https://www.spyder-ide.org/). Adjust **Input Folder and Files** section to refer to your local directory. Description of each variable:
> **PROJ_HOME** : main location of your run
>
> **SYS_APP** : location of your FnD3D folder, I suggest to always copy the folder under PROJ_HOME folder
>
> **D3D_HOME** : location of your Delft3D Flexible Mesh executables
>
> **gdal_loc** : location of your gdal installation, usually it is located inside your python virtual environment folder
>
> **JAVA_Exe** : location of your Java Emulator
>
> **Mangr_SHP** : mangrove shapefile name
>
> **D3D_Model** : folder name of your D3D run, always put under `Model-Execute\D3DFM`
>
> **D3D_Domain** : name of your DFM domain

Additionally, you can adjust another variable, following your specific model needs.

**EPSG_Project**, **x_res**, **y_res**, **LLWL**

## 2.2. Coupling Run
After conducting the Preparation Run, you should have the first MFON simulation and the outputs.

You can copy-paste similar files and folders location of your local run as in `01_prep_run.py` to `02_DFMFON_run.py`.

Complete the information of `config_xml` and `mdu_file`.

Congratulations, you will be ready for the full run of DFMFON model.



