# Drug parsers
The scripts that parse the Chembl drugs and compounds into our own RECOVER ids, and the scripts that align other sources of dtis to our own.

## Setup
Initially running `python create_initial_db_from_chembl.py` will convert ChEMBL drugs and compounds into our database of drugs. It will combine different molecules in chembl into one, using parsed smiles and salt hierarchies.

Afterwards the scripts from the dti folder can be run to map DTIs from Chembl, Drug repurposing hub and others to this reference:
* `python parse_chembl.py` - This parses the DTI as opposed to the actual compounds
* `python parse_drug_repurposing_hub.py`

## Molecule mapper
Molecule mapper is a package that is used to map molecules from a novel data sources to Recover ids.
