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
Numbers simply correspond to the number of carbon atoms in a linear alkyl chain. These can be augmented with SMILES a little, e.g. add an atom in place:

| P-***6***-0 | P-***2O2O***-0 | P-***2(C#C)3N2(C=C)3(C#C)2O***-0 |
|-------------|----------------|----------------------------------|
| ![download](https://github.com/user-attachments/assets/3dfac366-7148-47c8-b4a3-dcc45674e287) | ![download](https://github.com/user-attachments/assets/da35c8d0-1870-4a30-b829-49069f088f17) | ![download](https://github.com/user-attachments/assets/5d3cfacd-bc47-4bab-bd70-5897e4f8893f)

Letter-based fragments begin with an upper case letter, and include subsequent lower case letters up to the following upper case letter. These are used for both the core and side groups. 

So, a core wrriten as PPnN is read as "P" + "Pn" + "N", corresponding to 1,4-benzene, 2,5-benzonitrile, 2,5-pyridine. 

The core is given as the first part of the name, and the side groups are given in the second and third parts. It is prefrable to encode as much information in the core as possible, rather than using SMILES/SMARTS to augment the side groups. Some groups are used in both the core and side portions; they might have the same letter (e.g. ester / carboxylic acid as Z or Zl) or different (internal alkyne as "T" in core, terminal alkyne as "A" in side groups). 

You can generate these lists / images using the ```show_fragments``` function:

```cores, core_labels, core_img, ends, end_labels, end_img = show_fragments()```
### Core Units:
![download](https://github.com/user-attachments/assets/ad4166b1-a2ef-4086-9808-1ea35de30b95)
<br>
### Terminal Units:
![download](https://github.com/user-attachments/assets/d01c7f8d-371c-4129-bf85-b26477dd4d86)
<br>

## Tips and Tricks

You can insert SMILES strings  into the spacers to alter the structure; simple elements usually, but you can do more complex things (but that doesn't mean you should): <br>
| PP-6O-N | DUZGGQG-2S-N | PUMP-5Te-T | UUZG-4-0(C=CC#CC#N) |
|---------|--------------|------------|---------------------|
| ![download](https://github.com/user-attachments/assets/2ad23656-1c25-4ece-be45-651083e899e4) | ![download](https://github.com/user-attachments/assets/970d5080-5ee1-4724-aaa6-b4641bb12cf2)  | ![download](https://github.com/user-attachments/assets/5af47d69-e20b-4dd2-9346-911269df3957) | ![download](https://github.com/user-attachments/assets/35ca7e1a-af8b-441c-baa5-eb4a9e8c0b6b) |

<br>


There are different ways to name the same structure: <br>```smi = nts.process_msn('$STRING$')``` 
| VPPT-3-N | PPT-3(C=C)-N | PPV-0(N#CC#C)-3 | PP-3(C=C)-0(C#CC#N) |
|----------|--------------|-----------------|---------------------|
| ![download](https://github.com/user-attachments/assets/38d15b2f-b66c-4161-bb27-ae4bdaa97b08) | ![download](https://github.com/user-attachments/assets/5842a57f-79b3-4998-8252-d14720b80bc0) | ![download](https://github.com/user-attachments/assets/271e6fe2-c2fd-4e0b-acdd-147e9ae368b1) | ![download](https://github.com/user-attachments/assets/c3398dc1-c8e6-400c-8d16-08dcb224056d) | 
 
<br>
Generally, the best name is the shortest one that uses the least characters.

You can add new fragments to the JSON quite trivially; let's walk through how we'd add 3,6-substituted Tropolone:
*  We start with the SMILES string for Tropolone: ```C1=CC=C(C(=O)C=C1)O```
*  We add a connectivity marker for our "left side" attachment in the 3-position: ```C1=CC([X])=C(C(=O)C=C)1)O```
*  Analogously, we add a second connectivity marker for our "right side": ```C1=CC([X])=C(C(=O)C=C([Y])1)O```
*  We remove the digits for ring numbering and replace them with placeholders (& == 1, ^ == 2, % == 3 etc): ```C&=CC([X])=C(C(=O)C=C([Y])&)O```
*  We add to the JSON with an appropriate identifier (above): ```"Rop": "C&=CC([X])=C(C(=O)C=C([Y])&)O",```


## Questions
Q: Why use [X], [Y], &, ^, etc. in SMILES templates? These aren't SMILES friendly.<br>
A: These are placeholders that allow programmatic ring numbering and attachment logic; especially for nested, multi-ring systems. This lets us preserve the orientation of the rings too. You can "flip" the orientation of many rings by adding "l" after the fact.

Q: Can I extend this?
A: Yes! You can add new cores or terminal groups by modifying frags.json. 

Q: Does it handle stereochemistry?<br>
A: Basically, no. But you can pass this with [C@H] and [C@@H] in the chain/side parts, if you want.

Q: Why doesn't it have fragment XYZ?<br>
A: This is a tool to quickly decode things in patents, and to build datasets. If things are missing, they can be added to the JSON. 

Q: Why doesn't it support passing SMILES directly?
A: It could be made to do it, but why bother? By the time you've typed (or generated) the SMILES string for the fragment of interest you might as well have just done the entire molecule. As above, you can update the fragment list JSON file _easily_.

