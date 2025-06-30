# name_to_struct.py - convert "structure names" into smiles etc.
# why? Easy to curate a dataset if you can just type the molecules out

#import what we need
import numpy as np
from rdkit import Chem
import re
import json
import argparse
import sys

'''
Example usage:

in python
	load name_to_struct as nts:
	smi = nts.process_msn('DUZGU-3-F')
	print(smi)
	>> "CCCC0COC(c2cc(F)c(C(=O)Oc4cc(F)c(c6cc(F)c(F)c(F)c6)cc4)c(F)c2)OC0"

As a standalone:
	python name_to_struct.py "PP-5-N"
    >> PP-5-N → CCCCCc0ccc(c2ccc(C#N)cc2)cc0
    
On a list in a file:
    python name_to_struct.py mylist.txt
    >> PP-5-N → CCCCCc0ccc(c2ccc(C#N)cc2)cc0
    >> PPP-5-N → CCCCCc0ccc(c2ccc(c4ccc(C#N)cc4)cc2)cc0
    >> DUZGU-3-F → CCCC0COC(c2cc(F)c(C(=O)Oc4cc(F)c(c6cc(F)c(F)c(F)c6)cc4)c(F)c2)OC0
'''

# useful functions
def load_frags(filename = "frags.json"):
    # load the fragments from json and return as python dicts
    with open(filename, 'r') as f:
        frags = json.load(f)
    return frags['core_frags'], frags['end_frags']
    
def parse_string(string):
    return string.split('-')

def split_alpha_num(s: str) -> list[str]:
    """
    splits a string into runs of digits or letters.
    E.g. "3Se3" → ["3","Se","3"]. Handles brackets for branching
    Used for building smiles, later
    """
    pattern = r'\d+|[A-Za-z]+\([^\)]+\)|[A-Za-z]+|\([^\)]+\)'
    return re.findall(pattern, s)
    
def has_digits(s: str) -> bool:
    # just return true if string has a digit
    return any(char.isdigit() for char in s)

def split_by_caps(s: str) -> list[str]:
    # just split a string by caps, e.g. "PPlP" returns "P", "Pl", "P"
    # preserves brackets and their contents.
    return re.findall(r'[A-Z][a-z]*\(\d+\)?|[A-Z][a-z]*', s)

def format_ring_number(n: int) -> str:
    # make sure ring number obeys smiles syntax (<9, int; >=10, %int)
    return str(n) if n < 10 else f"%{n}"

def clean_smi_string(string):
    # function that will just clean a string:
    string = string.replace('((','(').replace('))',')') # remove double brackets
    string = string.replace('()','')
    return string
    
def update_core_numbers(processed_frag: str, count: int = 0) -> tuple[str, int]:
    """
    replace each placeholder symbol (&, ^, $, etc.) in matched pairs
    with the same ring number.
    """
    placeholders = ['&', '^', '$', '*', '~', '!', '@',]
    for ph in placeholders:
        ph_locs = [m.start() for m in re.finditer(re.escape(ph), processed_frag)]
        if len(ph_locs) % 2 != 0:
            raise ValueError(f"Unmatched ring placeholder: {ph} appears {len(ph_locs)} times")

        for i in range(0, len(ph_locs), 2):
            ring = format_ring_number(count)
            processed_frag = processed_frag.replace(ph, ring, 1)
            processed_frag = processed_frag.replace(ph, ring, 1)
            count += 1
    return processed_frag, count

def reverse_smiles_preserve_brackets(smi: str) -> str:
    #aim is to split a smiles string into tokens, preserving bracket order...
    # simply reversing the string (smi[::-1]) gives errors for things like C(=O)O >> O)O=(C!
    pattern = re.compile(r'[^\(\)]+(?:\([^\)]+\))*')
    
    tokens = pattern.findall(smi)

    if len(tokens) == 1: # no brackets, do it old way:
        reversed_smi = ''.join(tokens)[::-1]
        
    else: 
        tokens.reverse()
        reversed_smi = ''.join(tokens)
    
    return reversed_smi

def find_insertion_mode(processed_frag: str, patterns: list[str]) -> str:
    for pattern in patterns:
        if pattern in processed_frag:
            return pattern
    raise ValueError(f"No mode match found in: {processed_frag}")

def insert_string(string, insertion, index):
	# simply stick a string (string) inside another string (insertion) at the specified place (index)
    return string[:index] + insertion + string[index:]
	
# processing functions
def process_core(core_string, core_frags, smi = '', count = 0):
    core_string = split_by_caps(core_string)
    for n, a in enumerate(core_string):
        chain_str = None
        
        if '(' in a and ')' in a:
            # logic for adding a lateral THING
            a, chain_str = re.match(r'([A-Z][a-z]*?)\(([^()]*)\)', a).groups()
        
        processed_frag = (core_frags[a])

        # this loop handles the case ifwe have a lateral chain
        # we try to find a carbon atom on which to hang it
        if chain_str is not None:
            mode = find_insertion_mode(processed_frag, [
                'cc([Y])', 
                '([Y])cc', 
                '[Y]c', 
                'c[Y]', 
                'cc','n','c'
            ])
            lateral_chain, _ = process_chain(chain_str, smi='', mode=mode, count=0)
            processed_frag = insert_string(processed_frag, f"{lateral_chain}", processed_frag.find(mode))
    
        processed_frag, count = update_core_numbers(processed_frag, count)
        # depending on the position of the element of the name, its either:
        # core-chain/terminal-chain/terminal
        if n == 0:
            smi+=processed_frag
        elif n != 0:
            smi = smi.replace('[Y]',processed_frag.replace('[X]',''))
        if '&' in core_frags[a]:
            count+=1 # increment ring number counter
        if '^' in core_frags[a]:
            count+=1 # increment ring number counter
    return(smi, count)

def process_termini(term_string, end_frags, smi, mode = '[Y]', count = 0):
	# here, we'll compare our terminal string with the list of end_fragments to see what it is
    processed_frag = end_frags[term_string].replace('&',str(count))
    if mode == '[X]':
        processed_frag = reverse_smiles_preserve_brackets(processed_frag)
    smi = smi.replace(mode, processed_frag)
    if '&' in processed_frag:
        count+=1 
    return (smi, count)

def process_chain(chain_str: str, smi, mode = '[X]', count: int = 0) -> str:
    # process the terminal chain
    processed_frag = ""
    bits = split_alpha_num(chain_str)
    for b in bits:
        if b.isdigit():
            processed_frag += "C" * int(b)
        elif b.startswith("(") and b.endswith(")"):
            if mode == '[X]':
                b = b.replace('(','').replace(')','')
            processed_frag += b
        elif len(b) == 1 and b.isupper():
            processed_frag += b
        else:
            processed_frag += f"[{b}]"

    if mode not in ['[X]', '[Y]']:
        smi = insert_string(smi,"(" + processed_frag + ")", smi.find(mode))
    else:
        smi = smi.replace(mode, processed_frag)
    return smi, count # but we don't update count?
	
def process_msn(string):
    core_frags, end_frags = load_frags()
    parsed_string = parse_string(string)
    
    for n, ps in enumerate(parsed_string):
        mode = "[X]" if n == 1 else "[Y]" # set the mode incase the chain ins't in the expected position.
        if n == 0:
            smi, count = process_core(ps, core_frags, smi='', count=0)
     
        elif not has_digits(ps):
            # means its a terminus
            smi, count = process_termini(ps, end_frags, smi=smi, count=count, mode=mode)
        
        else:
            # means its a chain
            smi, count = process_chain(ps, smi=smi, mode = mode, count=count)
    
    smi = clean_smi_string(smi) # do a final cleanup
    return smi

# Show fragments function(s)
def show_fragments():
    core_frags, end_frags = load_frags()
    cores, core_labels = [], []
    ends, end_labels = [], []
    
    for key in core_frags:
        smi, _ = update_core_numbers(core_frags[key], count = 0)
        smi = smi.replace('[X]', '[*]')
        smi = smi.replace('[Y]', '[*]')
        smi = smi.replace('()','')
        
        core_labels.append(key)
        cores.append(Chem.MolFromSmiles(smi))
        
    core_img = Chem.Draw.MolsToGridImage(cores, 
                                         legends=core_labels, 
                                         maxMols = 100)
        
    for key in end_frags:
        smi = "[*]" + end_frags[key]
        smi = smi.replace('&','1')
        end_labels.append(key)
        ends.append(Chem.MolFromSmiles(smi))
        
    end_img = Chem.Draw.MolsToGridImage(ends, 
                                         legends=end_labels, 
                                         maxMols = 100)
        
    return cores, core_labels, core_img, ends, end_labels, end_img
cores, core_labels, core_img, ends, end_labels, end_img = show_fragments()

def main():
    parser = argparse.ArgumentParser(description="Convert structure names to SMILES strings.")
    parser.add_argument("input", help="Single structure name or path to .txt/.csv file with names")
    parser.add_argument("--out", help="Output file to write SMILES (optional)", default=None)
    
    args = parser.parse_args()

    input_val = args.input

    # if input is a file
    try:
        with open(input_val, 'r') as f:
            names = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        # if not a file — treat as single structure name
        names = [input_val]

    results = []
    for name in names:
        try:
            smi = process_msn(name)
            results.append((name, smi))
        except Exception as e:
            print(f"Failed to parse '{name}': {e}", file=sys.stderr)
            results.append((name, None))

    if args.out:
        with open(args.out, 'w') as f:
            for name, smi in results:
                f.write(f"{name},{smi or 'ERROR'}\n")
        print(f"Results written to {args.out}")
    else:
        for name, smi in results:
            print(f"{name} → {smi or 'ERROR'}")

if __name__ == "__main__":
    main()
