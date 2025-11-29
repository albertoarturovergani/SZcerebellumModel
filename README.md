# **SZcerebellumModel**

A computational framework for simulating **schizophrenia-like structural degeneration** in a biophysically detailed **cerebellar microcircuit**.

This repository extends the standard *cerebellum* + *BSB (Brain Scaffold Builder)* pipeline by introducing:

- Progressive neuronal **atrophy**  
  (*morphological shrinking, pruning, thinning*)
- Functional **resilience mechanisms**
- Automated pipelines for **building → simulating → analyzing** cerebellar networks
- Full support for:  
  **Mossy Fibers, Granule Cells, Golgi, Basket, Stellate, Purkinje Cells, DCN, IO**

---

## **Repository Structure**

This repository extends the standard `cerebellum/` framework by adding two main directories:

### **`script_sh/` — Shell Automation Scripts**
This folder contains all **bash pipelines** required for automated processing:

- morphological atrophy generation  
- network building (BSB)  
- simulation setup and execution  
- batch processing for multiple degeneration levels  
- post-simulation analysis orchestration  
- full end-to-end workflows (build → simulate → analyze)

---

### **`script_py/` — Python Processing Modules**
This folder includes all **Python tools** used throughout the workflow:

- morphology scaling, shrinking, pruning
- synaptic scaling  
- YAML generation for BSB network definitions  
- simulation launchers  
- detailed multi-level analyses (firing rate, TE, PLV, LFP, structure, etc.)  
- utilities for summarizing results across all degeneration levels

---


# **1. Prerequisites**

This repository **must be placed inside the `cerebellum/` directory**.  
Before using it, install the following components.

---

## **1.1 Brain Scaffold Builder (BSB)**

🔗 https://github.com/bsb-team/bsb

Typical installation:

```bash
pip install bsb-core
pip install bsb-cereb
```
## **1.2 cerebellum**

Official repository:  
🔗 https://github.com/bsb-team/cerebellum

---

## 2. Prepare folders

To enable automatic integration with the cerebellum pipeline,
copy or move the two directories script_sh and script_py into your local cerebellum/ installation:
```bash
cp -r script_sh   cerebellum/
cp -r script_py   cerebellum/
```
or
```bash
mv script_sh script_py cerebellum/
```
After this step, your folder structure should look like:

cerebellum/
    ├── cerebellum/
    ├── script_sh/
    ├── script_py/
    └── ...


## 3. Test Usage
Inside the cerebellum folder, run a full end-to-end test (degeneration → network build → simulation → analysis), execute:
```bash
source scripts_sh/run_test.sh
```
## Aknowledgemnts
This work is support by The Virtual Brain Twin Project, which has received funding from the European Union's Research and Innovation Program Horizon Europe under grant agreement No 101137289.

