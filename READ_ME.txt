All our calculations and analyses were performed using the operative system Windows (64 bits), with an open source version of Python: Anaconda 3, conda version : 24.11.3, conda-build version : 24.9.0, python version : 3.12.7.final.0, solver : libmamba (default)
       virtual packages : __archspec=1=haswell
                          __conda=24.11.3=0
                          __cuda=10.1=0
                          __win=0=0

Our analyses combined Flux Balance Analyses (performed with the open source software CobraPy) with thermodynamic calculations (performed with the open source software eQuilibrator). 

The inter-operatibility of these software inside the same environment depends on the combination of compatible versions. Finding by yourself the right combination(s) could be challenging. Therefore, to facilitate reproducibility of our analyses, we strongly recommend the reproduction, in your computer, of the Python environment we created for these analyses. 

To facilitate this process, we created the file named "cobra_and_equilibrator.yml". To reproduce in your computer the environment we created for our analyses, place the file named "cobra_and_equilibrator.yml" in an easy to find folder. Launch Anaconda prompt and navigate to the folder where you placed the file "cobra_and_equilibrator.yml". Then just type in the Anaconda prompt:

conda env create -f cobra_and_equilibrator.yml

Once the environment is created, you should activate it, typing in the Anaconda prompt:

conda activate cobra_and_equilibrator

Once the environment "cobra_and_equilibrator" is active, launch Jupyter Notebook from the Anaconda prompt, typing:

jupyter notebook

Once Jupyter Notebook is launched, navigate to the folder in your computer where you previously saved the following files:

custom_plot_functions.py
equilibrator_custom_functions.py
equilibrator_custom_functions_my.py
input C13 data.xlsx
input physiologic data during labeling.xlsx
input_data_bioreactor.csv
input_MW_values.csv
input_kinetic_parameters_database.csv
input_metabolite_ranges_default.csv
input_proteomics.xlsx
script_01_processing_proteomic_and_physiology_data.ipynb
script_02_FBA_stoichiometric_model_maker.ipynb
script_03_plotting_labeling_data.ipynb
script_04_FBA_metabolic_fluxes_labeling.ipynb
script_05_MDF_metabolic_fluxes_labeling.ipynb
script_06_FBA_pathway_maker.ipynb
script_07_MDF_pathways.ipynb
script_08_FBA_generating_file_for_ECM.ipynb
script_09_ECM_calculations.ipynb

The successful running of a given script depends on the existance of files created with a previous script. Therefore, you should execute the provided scripts by numerical order, i.e., first the script named "script_01_processing_proteomic_and_physiology_data.ipynb", then the script named "script_02_FBA_stoichiometric_model_maker.ipynb" and so on. 