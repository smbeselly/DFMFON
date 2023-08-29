# DFMFON
[![Last Commit](https://img.shields.io/github/last-commit/smbeselly/DFMFON/commits/main)](
https://github.com/smbeselly/DFMFON/commits/main)
[![GitHub release](https://img.shields.io/github/release/smbeselly/DFMFON)](https://GitHub.com/smbeselly/DFMFON/releases/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

DFMFON is an open-source software to mechanistically simulate the Mangrove and Hydromorphology development, written in Python.

This model is developed as part of the PhD research by Sebrian Beselly at [IHE Delft](https://www.un-ihe.org/department/coastal-and-urban-risk-resilience) and [Delft University of Technology](https://www.tudelft.nl/).

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


