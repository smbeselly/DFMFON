# DFMFON
[![Last Commit](https://img.shields.io/github/last-commit/smbeselly/DFMFON)](
https://github.com/smbeselly/DFMFON/commits/main)
[![GitHub release](https://img.shields.io/github/release/smbeselly/DFMFON)](https://GitHub.com/smbeselly/DFMFON/releases/)
This is a private repo of Delft3D-FM and MesoFON Coupling attempt. This is part of the PhD research by Sebrian Beselly in Monitoring and Modelling Mangrove-Mudflat Dynamics on a Prograding Delta

# What is DFMFON
In the face of climate change, sea level rise (SLR), increase precipitation certainly affects ecosytem and community therein.
Based on the report of IPCC in 2021, it is certain that the world has been experiencing an increase of temperature up to 1.5C.
That means, the effect of climate change is not negligible.
The rise of sea level is can not be ignored. Several reports has mentioned that several area of Earth has been experiencing SLR.
In response to that, we should investigate what kind of countermeasure that is sustainable for long run yet safe for all of the ecosystem.
Several research have shown how vegetation is dynamically adapt to sea level change. e.g., salt marsh and mangrove.
However, research that investigates the role of mangrove to sea level change is not yet understood.
In this research we attempt to understand how individual mangrove trees respond the environmental forces, e.g., wave force, tidal circulation, and sediment transport.

Here, we are trying to couple the state-of-the-art hydro-morphodynamic model Delft3D Flexible Mesh and an individual-based mangrove model MesoFON.
The output of this coupling project is the dynamic interactions up to the individual level of mangroves to environmental changes.

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


