# Reservoir

RECOVER coalition (Mila, Relation Therapeutics)

This Recover repository is based on research funded by (or in part by) the Bill & Melinda Gates Foundation. The findings and conclusions contained within are those of the authors and do not necessarily reflect positions or policies of the Bill & Melinda Gates Foundation.

This package stores the primary data acquisition scripts and the parsing pipelines for the main types of data in the RECOVER knowledge graph. It also provides an api for getting parsed data. This is **not** supposed to be used on the fly, it's just to get required data and save it in your dataset generation projects. This is because it is slow as it's constantly interacting with csv files

## Setup

-   Some files with the csv, tab, txt, obo, gaf and tsv are bigger than the 100MB limit and are being tracked with Git LFS. **Please install Git LFS** from here https://git-lfs.github.com/  **before** cloning the repository.
-   In order to be able to use the package you need to run `python setup.py develop` within the Reservoir directory.

## Structure

### data

Folder where the data is stored.:

-   raw: The data as it was extracted from the original data source
-   parsed: The data parsed from the raw format

### individual_pipelines

Here we store the scripts that are used to extract the data from their original source,  where we got the flat files from. Some examples:

-   Chembl: The data is extracted from their SQL database. We keep track of the SQL scripts used to extract this data
-   Drug repurposing hub: This data is in a flat file on their website. We just list the link to the file.

### parsers

Here we store the scripts used for parsing the raw data, and subsequently storing it in the parsed folder

## Api

This is the api for gaining access to parsed data without using the actual files

### Drugs and DTIs basics

    import reservoir as rsv

    # returns all the drugs in the db with some metadata
    drugs = rsv.get_drugs()

    # Get the available DTI files
    rsv.available_dtis()

    # Get the chembl binding dtis
    rsv.get_dtis('chembl_dtis.csv')

### Drugs advanced

    import reservoir as rsv

    # get the recover ids for some drugs using name, smiles or other db identifiers
    drugs_to_search = [
      {"drug_name": "Viagra"},
      {"input_smiles": "C[C@@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@@]3(F)[C@@H](O)C[C@]2(C)[C@@]1(O)C(=O)CO"},
      {"pubchem_cid": "5362132"},
      {"drug_name": "bla"}]
    mapped_drugs = rsv.map_drugs(pd.DataFrame(drugs_to_search))

    # get the dtis in chembl dtis for these drugs
    dtis = rsv.get_dtis("chembl_dtis.csv", recover_ids=set(mapped_drugs['recover_id']))

### Proteins and PPI basics

    import reservoir as rsv

    # returns all the proteins with inter db identifiers
    proteins = rsv.get_proteins()

    # get the available ppi files
    rsv.available_ppis()

    # get the huri ppis
    rsv.get_ppis("huri.csv")

### Proteins advanced

    import reservoir as rsv

    # get ppi involving the ACE2 and NFE2L2 proteins
    rsv.get_ppis("biogrid.csv", gene_hgnc_ids=['ACE2', 'NFE2L2'])

### Gene ontology - basics

    import reservoir as rsv

    # get all annotations
    rsv.get_gene_ontology()

    # annotations can be filtered by a few criteria
    rsv.get_gene_ontology(gene_hgnc_ids=['ACE2'], aspect='biological_process', max_distance_from_root=2, evidence_category='experimental evidence')

    # the types of evidence categories available
    rsv.get_gene_ontology_evidence_categories()

### Gene ontology - embeddings

You can create simple embeddings for proteins using Gene Ontology. The `rsv.get_gene_onotology_embeddings()` function does this by
identifying all of the unique annotations (N) that exist at a certain level in the hirearchy. Each protein is then represented by a length N vector of 1s and 0s indicating if a certain annotation is present or not

    # get embeddings at level 2 in the tree using annotations that have experimental backing
    rsv.get_gene_onotology_embeddings(distance_from_root=2, evidence_category='experimental evidence')

### Covid
The Covid proteins can be accessed using:
```python
rsv.get_covid_proteins()
```

The viral-human PPIs can be accessed using:
```python
rsv.get_ppis("covid_krogan_ppi.csv")
rsv.get_ppis("covid_liang_ppi.csv")
```
Specific information about the PPIs can be found here https://github.com/RECOVERcoalition/Reservoir/tree/master/reservoir/data/parsed/ppi

### Drug combos from https://drugcomb.fimm.fi/
Drug combos are organised by blocks, and each block has associated metadata like cell line used, time, the original study, the number of measurements for the combos and for the monotherapies

The mono-therapies are distinguished from each other using the row col names as per their original source.

First, let's say you want to get all of the blocks that have 9 measurements for the combo and 4 for the row and col monotherapies
```python
import reservoir as rsv
blocks = rsv.get_specific_drug_combo_blocks(combo_measurements=9, mono_row_measurements=4, mono_col_measurements=4)
```

Now you want to get the mono-therapy data for those blocks
```python
mono_data = rsv.get_drug_combo_data_monos(block_ids=blocks['block_id'])

# Show the data for the row drug
print(mono_data[['block_id', 'conc_r', 'inhibition_r']])

# Show the data for the col drug
print(mono_data[['block_id', 'conc_c', 'inhibition_c']])

# These can be stacked as they are numpy arrays
import numpy as np
np.stack(mono_data['inhibition_r'])
```

Now the combo data
```python
combo_data = rsv.get_drug_combo_data_combos(block_ids=blocks['block_id'])

# The inhibition measurements from the combination matrix are here as a list
combo_data[['block_id', 'inhibitions']]

# The concentration for the combination matrix are here, each tuple is  (row drug, column drug). They map by index to the inhibitions
combo_data[['block_id', 'concentration_pairs']]

# Synergy scores
combo_data[['block_id', 'css_ri', 'synergy_zip', 'synergy_bliss', 'synergy_loewe', 'synergy_hsa', 'S']]
```

Use argument version=1.5 or 1.4 to switch between new and old versions of DrugComb
```python
# DrugComb 1.5 combos with medium QC filters
combo_data = rsv.get_drug_combo_data_combos(qc_filtering='medium')

# Almanac blocks DrugComb 1.4
almanac = rsv.get_specific_drug_combo_blocks(study_name='ALMANAC')
```
