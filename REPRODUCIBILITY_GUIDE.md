# Reproducing Our Computational Analyses

## System and Environment

All calculations and analyses were performed on **Windows (64-bit)** using **Anaconda 3** with the following specifications:

- **Conda version**: 24.11.3  
- **Conda-build version**: 24.9.0  
- **Python version**: 3.12.7.final.0  
- **Solver**: libmamba (default)  

**Virtual packages:**
```
__archspec=1=haswell
__conda=24.11.3=0
__cuda=10.1=0
__win=0=0
```

We used the open-source packages **CobraPy** (for Flux Balance Analysis) and **eQuilibrator** (for thermodynamic calculations). Compatibility between these tools requires specific package versions, so we strongly recommend recreating our Python environment.

---

## Required Files

Create a folder on your computer and place the following files inside:

```
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
```

---

## Option 1: Using the Predefined Environment File (Recommended)

To ensure all software dependencies are correctly installed, use the environment file we created.

### Steps:

1. **Download** the file `cobra_and_equilibrator.yml` and place it in an easy-to-find folder.
2. **Open** the **Anaconda Prompt**.
3. **Navigate** to the folder where the `.yml` file is located. For example:
   ```
   cd path\to\your\folder
   ```
4. **Create the environment** by running:
   ```
   conda env create -f cobra_and_equilibrator.yml
   ```
5. **Activate** the new environment:
   ```
   conda activate cobra_and_equilibrator
   ```
6. **Launch** Jupyter Notebook:
   ```
   jupyter notebook
   ```

7. In the browser window that opens, navigate to the folder containing the files and run the scripts **in numerical order**:

   - `script_01_processing_proteomic_and_physiology_data.ipynb`
   - `script_02_FBA_stoichiometric_model_maker.ipynb`
   - ...
   - `script_09_ECM_calculations.ipynb`

> ⚠️ Note: Some scripts generate input files for the next ones. Skipping or reordering them may result in missing dependencies or errors.

You can find the necessary input files and scripts at our GitHub repository:  
👉 [https://github.com/kolavarria/NTS_pathway](https://github.com/kolavarria/NTS_pathway)

---

## Option 2: Manually Creating the Environment

If you prefer to set up the environment yourself, install the following packages **in this order**:

```bash
conda create -n cobra_and_equilibrator python=3.9
conda activate cobra_and_equilibrator

pip install equilibrator-api==0.4.5.1 equilibrator-cache==0.4.3 equilibrator-pathway==0.4.4

conda install -c conda-forge cobra
pip install notebook==6.5.6

python -m ipykernel install --user --name=cobra_and_equilibrator --display-name "Python (cobra_and_equilibrator)"

conda install pandas=1.5.3 sqlalchemy=1.4.49
conda install -c conda-forge plotly=5.15.0
pip install -U kaleido
pip install pymupdf
```

> ⚠️ Caution: We have not tested other versions of these packages. Proceed at your own risk if using alternatives.
