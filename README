# Revarcine: a reverse vaccinology pipeline

This project was created and is maintained by Omixlab team.

---

### **About**

This module characterizes the user through Signal Peptide, Signal Peptide Cut Position, Subcellular Location, Beta Barrel and Alpha Helix prediction. In this context, **Signal Peptide** means predict the presence or absence of the signal peptide structure in the amino acid sequence. The **Signal Peptide Cut Position** is a model that predicts the amino acid that represents the cut position of the signal peptiode in sequences that have this structure. **Subcellular Location** is used to determine the subcellular location of the protein. The **Beta Barrel** analysis indicates the presence or absence of this structure. The **Alpha Helix** predict the presence or absence of this structure.

### Project structure

```
└───revarcine
    ├───data/                  ---> data input for revarcine prediction
    │   ├───input/             ---> input folder
    │   └───output/            ---> predict output folder
    └───exe/                   ---> execution structure files
    │   ├───model/             ---> h5 model file
    │   └───tokenizer/         ---> pikle model file
    ├───environment.yml        ---> conda environment file
    ├───LICENSE                ---> license for code distribution
    ├───Makefile               ---> Makefile to centralize the main analysis steps in "rules"
    ├───README.md              ---> Project README
    └───requirements.txt       ---> python requirements
```

#### Parameters
| Param  | Obligatory | Description  |	Example	|
| :-- | :--: | :-- | :-- |
| **INPUT_NAME** |&#10003;| Input file's path .fasta file |	<code>micobacterium.fasta</code>	|
| **INPUT_GRAM** | &#10003; |  Specify how gram classification the bacteria genome is G+ ou G- | <code>G+</code> |
| **ARGPRED** | &#10003; |  Specify which analysis will be run | <code> 10010 </code> |

---

### **Analysis Combination (**ARGPRED**)**
Each analysis correspond to a binary flag in a list which its index is linked to the model's name.

**Models sequency:** **signal peptide, signal peptide cut position, subcellular location, beta barrel, alpha helix**

**E.g**:

1. **Default - 11111** : Execute all analysis wich includes: signal peptide, subcellular location, signal peptide cut position, beta abrrel and alpha helix .

3. **10000** : to run just Sigal Peptide analysis.

4. **00100** : to run just Subcellular Location analysis.

5. **00011** : to run Beta Barrel and Alpha Helix analysis.

---
### **Base Input for all analysis**

Please provide a base with the .fasta format, use the exemple bellow:

>O50499
MTIRAGSLDRRTLLRGAIATAAMGSFAVACSSPSSEDKESDSGPKGEKSANNPFGAAANSTVEAAIFDGGYGTDYVDYANQVLGSQVKGLKVQVKPVVDIAPQLQPRFVGGNPPDLIDNSGEDQIGFLGILDQLEELDDLFEASTYEGKKIADIVYPGVKDPGTFKDKFVALNYVMTVYGVWYSKTLFEENGWTPPKTWDEALDLGQEAKKKGKYLFVHGKEAATYYRTLLIDSAIKEGGDEVRLALENLEKGCWSHPAVQG

---

### Install

###### 1. Open a shell and clone this repo 

```bash
git clone git@github.com:omixlab/revarcine.git
```

###### 2. Go to **revarcine/** folder.

###### 3. Create revarcine environment.

###### 4. Install make.
* **For windows:**
1- Open PowerShell (Windows + x)
2- Install Chocolatey

```bash
Set-ExecutionPolicy Bypass -Scope Process -Force; iex ((New-Object System.Net.WebClient).DownloadString('https://community.chocolatey.org/install.ps1'))
```

3- Verify Chocolatey install

```bash
choco
```

4- Make install

```bash
choco install make
```

* **For linux:**

```bash
sudo apt install make
```

###### 5. Conda revarcine environment activation.

```bash
conda activate revarcine
```

###### 6. Conda revarcine environment activation.

```
make upgrade
```

**Remember to activate the environment every time you run the revarcine**

---
### Usage

### 1. Run
```python
python-W ignore revarcine.py INPUT_NAME 'INPUT_GRAM' 'ARGPRED'

```
### Example:
```python
python -W ignore revarcine.py exemple.fasta 'G+' '11000'

```
