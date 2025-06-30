# name_to_struct

A lightweight Python utility to convert "structure names" for liquid crystals, basically — a shorthand chemical naming system — into valid SMILES strings. The use case for this is quickly curating molecular datasets by enabling me/us to type concise names rather than draw full chemical structures, or write valid SMILES strings.

---

## Features

- Converts molecule shorthand (e.g. `PP-5-N`) into canonical SMILES strings. This shorthand is widely used in the patent literature.
- We break down a structure into three parts: "core" - "left side" - "right side". The core is made up of rings etc, (e.g. P = phenyl); the sides can chains (e.g. 5 = pentyl) or functional groups (e.g. N = nitrile). 
- Handles:
  - Arbitrary ring systems with placeholder numbering.
  - Branched/lateral substituents via embedded parentheses (e.g. P(2) is ethyl benzene; M(1) is methylpyrimidine and so on.
  - Chain lengths can be complex by using brackets to indicate branching and passing SMILES strings, e.g. 2(C)1(F)(Cl) gives a 1-chloro-1-fluoro-2-methylpropyl chain.
  - Many widely used terminal groups.
  - Easy to include new structures via included JSON file.
- Works on single names or batches of names from a file. You can run it via the commandline, or import and use in a notebook.
- Designed to integrate with RDKit pipelines.

---

## Requirements

- Python 3.7+
- RDKit
- NumPy

## Usage
### Jupyter Notebook:
```
import name_to_struct as nts
from rdkit import Chem

smi = nts.process_msn('DUZGU-3-F')
mol = Chem.MolFromSmiles(smi)
mol
```
![download](https://github.com/user-attachments/assets/7caf31a5-6ac2-49e1-aa51-2fbb665fecc3)
### Terminal 
On a single molecule:
```python name_to_struct.py "DUZPUQU-3-F"```
![image](https://github.com/user-attachments/assets/f2ddc411-49f7-4c50-80bc-16887d66f6e9)

On many molecules (put the names in a text or csv file; one entry per line):
```python name_to_struct.py mylist.txt```
![image](https://github.com/user-attachments/assets/b36db25d-2aaf-4a25-8a6b-255e32a6f851)

## Coded Fragments:
You can generate these lists / images using the ```show_fragments``` function:

```cores, core_labels, core_img, ends, end_labels, end_img = show_fragments()```
### Core Units:
![download](https://github.com/user-attachments/assets/76594044-dc20-43bd-8e1a-587a9e915aaf)<br>
### Terminal Units:
![download](https://github.com/user-attachments/assets/c7ea1515-9caf-49e0-9160-ca23a6af3ecf)<br>

## Questions
Q: Why use [X], [Y], &, ^, etc. in SMILES templates? These aren't SMILES friendly.
A: These are placeholders that allow programmatic ring numbering and attachment logic — especially for nested, multi-ring systems. This lets us preserve the orientation of the rings too. You can "flip" the orientation of many rings by adding "l" after the fact.

Q: Can I extend this?
Yes! You can add new cores or terminal groups by modifying frags.json. 

Q: Does it handle stereochemistry?
Basically, no. But you can pass this with [C@H] and [C@@H] in the chain/side parts, if you want.

