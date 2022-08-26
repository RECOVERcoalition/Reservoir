import pandas as pd
import reservoir as rsv
import numpy as np
import json
import urllib.request
import os.path
from multiprocessing import Pool


def get_data_api(block_id):
    """ Get assay data from drugcomb API
    """
    url = 'https://api.drugcomb.org/response/' + str(block_id)
    with urllib.request.urlopen(url) as response:
        html = response.read()
    encoding = response.headers.get_content_charset('utf-8')
    html_text = html.decode(encoding)
    json_text = json.loads(html_text)
    block = pd.DataFrame(json_text)
    return block


def get_concentration_pairs(rows):
    """ Function for extracting concentration pairs for each block
    """
    pairs = list(rows.apply(lambda row: tuple([row["conc_r"], row["conc_c"]]), axis=1))
    inhibitions = list(rows["inhibition"])

    return pd.DataFrame(
        [
            {
                "drug_row_RESERVOIR_id": rows.iloc[0]["drug_row_RESERVOIR_id"],
                "drug_row_smiles": rows.iloc[0]["drug_row_smiles"],
                "drug_col_RESERVOIR_id": rows.iloc[0]["drug_col_RESERVOIR_id"],
                "drug_col_smiles": rows.iloc[0]["drug_col_smiles"],
                "cell_line_name": rows.iloc[0]["cell_line_name"],
                "pairs": pairs,
                "inhibitions": inhibitions,
            }
        ]
    )


chunksize = 10000 # For multiprocessing

summary_data = pd.read_csv(rsv.RESERVOIR_DATA_FOLDER + "/raw/drug-combos/v1.5/summary_v_1_5.csv")

fname = rsv.RESERVOIR_DATA_FOLDER + "/raw/drug-combos/v1.5/blocks_data.csv"

if not os.path.isfile(fname):
    with Pool(100) as p:
        summary_data['block_data'] = p.map(get_data_api, summary_data['block_id'], chunksize)
    blocks = pd.concat(summary_data['block_data'].values)
    blocks.to_csv(fname, index=False)
    summary_data = summary_data.drop(columns='block_data')
else:
    blocks = pd.read_csv(fname)


# Get all unique drugs
row_drugs = (
    summary_data[["drug_row"]]
    .drop_duplicates()
    .rename(columns={"drug_row": "drug_name"})
)
col_drugs = (
    summary_data[["drug_col"]]
    .drop_duplicates()
    .rename(columns={"drug_col": "drug_name"})
)
unique_drugs = pd.concat([row_drugs, col_drugs]).drop_duplicates()
unique_drugs = unique_drugs.loc[~unique_drugs["drug_name"].isna()]

# Add drug information
# drugs file from https://api.drugcomb.org/drugs
drug_json = open(rsv.RESERVOIR_DATA_FOLDER + "/raw/drug-combos/v1.5/drugs")
drug_json_data = json.load(drug_json)
drug_info = pd.DataFrame(drug_json_data)
drug_info = drug_info.rename(columns={'dname': 'drug_name',
                                      'smiles': 'inputs_smiles',
                                      'cid': 'pubchem_cid'})

unique_drugs = unique_drugs.merge(drug_info, on=['drug_name'])
unique_drugs = unique_drugs.drop_duplicates(subset=['drug_name'])

# Map drugs to RESERVOIR IDs
unique_drugs['pubchem_cid'] = unique_drugs['pubchem_cid'].astype(object) # map_drugs fails if pubchem cid is float
mapped_drugs = rsv.map_drugs(unique_drugs)
mapped_drugs = mapped_drugs.loc[~pd.isna(mapped_drugs["RESERVOIR_id"])]
mapped_drug_data = rsv.get_drugs(filter=set(mapped_drugs["RESERVOIR_id"]))
mapped_drugs = mapped_drugs.merge(mapped_drug_data[["RESERVOIR_id", "smiles"]])

# Drug name to lower case to allow merging with mapped drugs
summary_data['drug_col'] = summary_data['drug_col'].str.lower()
summary_data['drug_row'] = summary_data['drug_row'].str.lower()

# Connect to dose data
drug_combos = summary_data.merge(
    mapped_drugs[["RESERVOIR_id", "smiles", "drug_name"]].rename(
        columns={"RESERVOIR_id": "drug_row_RESERVOIR_id", "smiles": "drug_row_smiles"}
    ),
    left_on="drug_row",
    right_on="drug_name",
)

drug_combos = drug_combos.merge(
    mapped_drugs[["RESERVOIR_id", "smiles", "drug_name"]].rename(
        columns={"RESERVOIR_id": "drug_col_RESERVOIR_id", "smiles": "drug_col_smiles"}
    ),
    left_on="drug_col",
    right_on="drug_name",
)

# Add block info to drug combos
drug_combos = drug_combos.rename(columns={'synergy_zip': 'synergy_zip_summary',
                                          'synergy_loewe': 'synergy_loewe_summary',
                                          'synergy_hsa': 'synergy_hsa_summary',
                                          'synergy_bliss': 'synergy_bliss_summary'})

drug_combos = drug_combos.merge(blocks, on='block_id')

# Extract monotherapy for the row side
mono_therapy_row = drug_combos.loc[
    drug_combos["conc_c"] == 0.0,
    [
        "block_id",
        "drug_row_RESERVOIR_id",
        "drug_row_smiles",
        "cell_line_name",
        "conc_r",
        "inhibition",
    ],
]
mono_therapy_row_by_block = (
    mono_therapy_row.groupby("block_id")
    .agg(
        {
            "drug_row_RESERVOIR_id": "first",
            "drug_row_smiles": "first",
            "cell_line_name": "first",
            "conc_r": list,
            "inhibition": list,
        }
    )
    .reset_index()
)
mono_therapy_row_by_block = mono_therapy_row_by_block.rename(
    columns={"inhibition": "inhibition_r"}
)
mono_therapy_row_by_block["conc_r"] = mono_therapy_row_by_block["conc_r"].apply(
    np.array
)
mono_therapy_row_by_block["inhibition_r"] = mono_therapy_row_by_block[
    "inhibition_r"
].apply(np.array)

# Extract monotherapy for the col side
mono_therapy_col = drug_combos.loc[
    drug_combos["conc_r"] == 0.0,
    [
        "block_id",
        "drug_col_RESERVOIR_id",
        "drug_col_smiles",
        "cell_line_name",
        "conc_c",
        "inhibition",
    ],
]
mono_therapy_col_by_block = (
    mono_therapy_col.groupby("block_id")
    .agg(
        {
            "drug_col_RESERVOIR_id": "first",
            "drug_col_smiles": "first",
            "cell_line_name": "first",
            "conc_c": list,
            "inhibition": list,
        }
    )
    .reset_index()
)
mono_therapy_col_by_block = mono_therapy_col_by_block.rename(
    columns={"inhibition": "inhibition_c"}
)
mono_therapy_col_by_block["conc_c"] = mono_therapy_col_by_block["conc_c"].apply(
    np.array
)
mono_therapy_col_by_block["inhibition_c"] = mono_therapy_col_by_block[
    "inhibition_c"
].apply(np.array)

# combine row and col monotherapy data
mono_therapy_combined = mono_therapy_row_by_block.merge(mono_therapy_col_by_block)
mono_therapy_combined = mono_therapy_combined[
    [
        "block_id",
        "cell_line_name",
        "drug_row_RESERVOIR_id",
        "drug_row_smiles",
        "conc_r",
        "inhibition_r",
        "drug_col_RESERVOIR_id",
        "drug_col_smiles",
        "conc_c",
        "inhibition_c",
    ]
]


# Extract combos and groupby block
combo_therapy = drug_combos.loc[
    (drug_combos["conc_c"] != 0.0) & (drug_combos["conc_r"] != 0.0)
]
combo_therapy_by_block = combo_therapy.groupby("block_id").apply(
    get_concentration_pairs
)
combo_therapy_by_block = combo_therapy_by_block.reset_index()

# Connect with synergy scores
combo_therapy_by_block = combo_therapy_by_block.merge(
    summary_data[
        [
            "block_id",
            "ic50_row",
            "ic50_col",
            "ri_row",
            "ri_col",
            "css_row",
            "css_col",
            "css_ri",
            "S_sum", "S_mean", "S_max",
            "synergy_zip",
            "synergy_bliss",
            "synergy_loewe",
            "synergy_hsa",
            "drug_row_clinical_phase",
            "drug_col_clinical_phase",
            "drug_row_target_name",
            "drug_col_target_name"
        ]
    ]
)

combo_therapy_by_block = combo_therapy_by_block[
    [
        "block_id",
        "cell_line_name",
        "drug_row_RESERVOIR_id",
        "drug_row_smiles",
        "drug_col_RESERVOIR_id",
        "drug_col_smiles",
        "pairs",
        "inhibitions",
        "ic50_row",
        "ic50_col",
        "ri_row",
        "ri_col",
        "css_row",
        "css_col",
        "css_ri",
        "S_sum", "S_mean", "S_max",
        "synergy_zip",
        "synergy_bliss",
        "synergy_loewe",
        "synergy_hsa",
        "drug_row_clinical_phase",
        "drug_col_clinical_phase",
    ]
]

combo_therapy_by_block = combo_therapy_by_block.rename(
    columns={"pairs": "concentration_pairs"}
)

# Add min and max synergy scores
min_synergy = drug_combos[["block_id",
                           "synergy_zip", 
                           "synergy_bliss",
                           "synergy_loewe",
                           "synergy_hsa"]].groupby('block_id').agg('min').reset_index()

new_col_names = list()
for col in min_synergy.columns:
    if col != 'block_id':
        col = col + '_min'
    new_col_names.append(col)

min_synergy.columns = new_col_names

max_synergy = drug_combos[["block_id",
                           "synergy_zip",
                           "synergy_bliss",
                           "synergy_loewe",
                           "synergy_hsa"]].groupby('block_id').agg('max').reset_index()

new_col_names = list()
for col in max_synergy.columns:
    if col != 'block_id':
        col = col + '_max'
    new_col_names.append(col)
    
max_synergy.columns = new_col_names

combo_therapy_by_block = combo_therapy_by_block.merge(min_synergy, on='block_id')
combo_therapy_by_block = combo_therapy_by_block.merge(max_synergy, on='block_id')

# Some blocks contain the inhibitions matrix twice. We need to clean that up.
doubled = combo_therapy_by_block.concentration_pairs.apply(
    lambda cp: len(set(tuple(p) for p in cp)) != len(cp) )
    
combo_therapy_by_block.loc[doubled,'concentration_pairs'] = combo_therapy_by_block.loc[
    doubled,'concentration_pairs'].apply(lambda c: c[::2])
combo_therapy_by_block.loc[doubled,'inhibitions'] = combo_therapy_by_block.loc[
    doubled,'inhibitions'].apply(lambda c: c[::2])


#%%export

combo_therapy_by_block.to_json(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/drug_combos/v1.5/combos.json"
)

mono_therapy_combined.to_json(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/drug_combos/v1.5/mono_therapy.json"
)
#%%

# add number of measurements from mono and combo therapies
mono_therapy_combined["mono_row_measurements"] = mono_therapy_combined[
    "inhibition_r"
].apply(len)
mono_therapy_combined["mono_col_measurements"] = mono_therapy_combined[
    "inhibition_c"
].apply(len)
combo_therapy_by_block["combo_measurements"] = combo_therapy_by_block[
    "inhibitions"
].apply(len)

# clean summary data
summary_data = summary_data[
        ["block_id", "study_name"]
]

summary_data = summary_data.merge(
    combo_therapy_by_block[["block_id", "cell_line_name"]]
)

summary_data = summary_data.merge(
    mono_therapy_combined[
        ["block_id", "mono_row_measurements", "mono_col_measurements"]
    ]
)
summary_data = summary_data.merge(
    combo_therapy_by_block[["block_id", "combo_measurements"]]
)

#%%export

summary_data.to_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/drug_combos/v1.5/summary_data.csv", 
    index=False
)

#%%Quality control

#filter out failed/low quality experiments
with Pool(100) as p:
    combo_therapy_by_block['inhibi_mean'] = p.map(np.mean, combo_therapy_by_block.inhibitions)

with Pool(100) as p:
    combo_therapy_by_block['inhibi_std'] = p.map(np.std, combo_therapy_by_block.inhibitions)

block_mask = \
    ( combo_therapy_by_block['inhibi_std'] > 0.05 ) & \
    ( combo_therapy_by_block["combo_measurements"] >= 9 ) & \
    ( combo_therapy_by_block['inhibi_mean'] >= -15 ) & \
    ( combo_therapy_by_block['inhibi_mean'] < 100 )

three_by_three_mask = list()
for c in combo_therapy_by_block['concentration_pairs']:
    col = pd.DataFrame(c)[0].drop_duplicates()
    row = pd.DataFrame(c)[1].drop_duplicates()
    if len(col) > 2 and len(row) > 2:
        three_by_three_mask.append(True)
    else:
        three_by_three_mask.append(False)

three_by_three_mask = pd.Series(three_by_three_mask)

block_mask = block_mask & three_by_three_mask
block_mask.index = combo_therapy_by_block.block_id

block_mask = block_mask.reset_index()
block_mask.columns = ['block_id', 'mask']

#%%export
block_mask[['block_id', 'mask']].to_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/drug_combos/v1.5/block_mask.csv",
    index=False
    )

#%%filter again with stricter thresholds
block_mask = \
    ( combo_therapy_by_block['inhibi_std'] > 0.05 ) & \
    ( combo_therapy_by_block["combo_measurements"] >= 9 ) & \
    ( combo_therapy_by_block['inhibi_mean'] >= -5 ) & \
    ( combo_therapy_by_block['inhibi_mean'] < 95 )


three_by_three_mask = list()
for c in combo_therapy_by_block['concentration_pairs']:
    col = pd.DataFrame(c)[0].drop_duplicates()
    row = pd.DataFrame(c)[1].drop_duplicates()
    if len(col) > 2 and len(row) > 2:
        three_by_three_mask.append(True)
    else:
        three_by_three_mask.append(False)

three_by_three_mask = pd.Series(three_by_three_mask)
block_mask = block_mask & three_by_three_mask

block_mask.index = combo_therapy_by_block.block_id

block_mask = block_mask.reset_index()
block_mask.columns = ['block_id', 'mask']

#%%export2
block_mask[['block_id', 'mask']].to_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/drug_combos/v1.5/block_mask_hq.csv",
    index=False
    )
