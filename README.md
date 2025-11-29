# **SZcerebellumModel**

A computational framework for simulating **schizophrenia-like structural degeneration** in a biophysically detailed **cerebellar microcircuit**.

This repository extends the standard *cerebellum* + *BSB (Brain Scaffold Builder)* pipeline by introducing:

- Progressive neuronal **atrophy**  
  (*morphological shrinking, pruning, geodesic reduction*)
- Synaptic **rewiring** driven by morphological degeneration  
- Functional **resilience mechanisms**
- Automated pipelines for **building → simulating → analyzing** cerebellar networks
- Full support for:  
  **Mossy Fibers, Granule Cells, Golgi, Basket, Stellate, Purkinje Cells, DCN, IO**

---

## **Repository Structure**

*(You may fill this section with folder descriptions, workflow diagrams, etc.)*

cartelle script_sh e script_py 

# **1. Prerequisites**

This repository **must be placed inside the `cerebellum/` directory**.  
Before using it, install the following components.

---

## **1.1 cerebellum**

Official repository:  
🔗 https://github.com/bsb-team/cerebellum

---

## **1.2 Brain Scaffold Builder (BSB)**

🔗 https://github.com/bsb-team/bsb

Typical installation:

```bash
pip install bsb-core
pip install bsb-cereb
```
## 2. Prepare folders

To enable automatic integration with the cerebellum pipeline,
copy or move the two directories script_sh and script_py into your local cerebellum/ installation:

cp -r script_sh   cerebellum/
cp -r script_py   cerebellum/

or

mv script_sh script_py cerebellum/

After this step, your folder structure should look like:

cerebellum/
    ├── cerebellum/
    ├── examples/
    ├── script_sh/
    ├── script_py/
    └── ...


## 3. Test Usage
Inside the cerebellum folder, run a full end-to-end test
(degeneration → network build → simulation → analysis), execute:

source scripts_sh/run_test.sh
