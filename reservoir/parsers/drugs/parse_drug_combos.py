import pandas as pd
import reservoir as rsv
import itertools

summary_data = pd.read_csv(rsv.RESERVOIR_DATA_FOLDER + "/raw/drug_combos/assay_v5.csv")
drug_combos = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/raw/drug_combos/drugcomb_data_v1.4.csv"
)
combo_scores = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/raw/drug_combos/DataTable5fb4ea7b3b1cb78467.csv"
)
drug_combos = drug_combos.astype({"drug_row_cid": "str", "drug_col_cid": "str"})

# get all unique drugs
row_drugs = (
    drug_combos[["drug_row", "drug_row_cid"]]
    .drop_duplicates()
    .rename(columns={"drug_row": "drug_name", "drug_row_cid": "pubchem_cid"})
)
col_drugs = (
    drug_combos[["drug_col", "drug_col_cid"]]
    .drop_duplicates()
    .rename(columns={"drug_col": "drug_name", "drug_col_cid": "pubchem_cid"})
)
unique_drugs = pd.concat([row_drugs, col_drugs]).drop_duplicates()
unique_drugs = unique_drugs.loc[~unique_drugs["drug_name"].isna()]

# map drugs
mapped_drugs = rsv.map_drugs(unique_drugs)
mapped_drugs = mapped_drugs.loc[~pd.isna(mapped_drugs["recover_id"])]
mapped_drug_data = rsv.get_drugs(filter=set(mapped_drugs["recover_id"]))
mapped_drugs = mapped_drugs.merge(mapped_drug_data[["recover_id", "smiles"]])

# connect to dose data
drug_combos = drug_combos.merge(
    mapped_drugs[["recover_id", "smiles", "pubchem_cid"]].rename(
        columns={"recover_id": "drug_row_recover_id", "smiles": "drug_row_smiles"}
    ),
    left_on="drug_row_cid",
    right_on="pubchem_cid",
)
drug_combos = drug_combos.merge(
    mapped_drugs[["recover_id", "smiles", "pubchem_cid"]].rename(
        columns={"recover_id": "drug_col_recover_id", "smiles": "drug_col_smiles"}
    ),
    left_on="drug_col_cid",
    right_on="pubchem_cid",
)

# extract monotherapy for the row side
mono_therapy_row = drug_combos.loc[
    drug_combos["conc_c"] == 0.0,
    [
        "block_id",
        "drug_row_recover_id",
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
            "drug_row_recover_id": "first",
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

# extract monotherapy for the col side
mono_therapy_col = drug_combos.loc[
    drug_combos["conc_r"] == 0.0,
    [
        "block_id",
        "drug_col_recover_id",
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
            "drug_col_recover_id": "first",
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
        "drug_row_recover_id",
        "drug_row_smiles",
        "conc_r",
        "inhibition_r",
        "drug_col_recover_id",
        "drug_col_smiles",
        "conc_c",
        "inhibition_c",
    ]
]
mono_therapy_combined.to_json(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/drug_combos/mono_therapy.json"
)


def get_concentration_pairs(rows):
    """ Function for extracting concentration pairs for each block
    """
    pairs = list(rows.apply(lambda row: tuple([row["conc_r"], row["conc_c"]]), axis=1))
    inhibitions = list(rows["inhibition"])

    return pd.DataFrame(
        [
            {
                "drug_row_recover_id": rows.iloc[0]["drug_row_recover_id"],
                "drug_row_smiles": rows.iloc[0]["drug_row_smiles"],
                "drug_col_recover_id": rows.iloc[0]["drug_col_recover_id"],
                "drug_col_smiles": rows.iloc[0]["drug_col_smiles"],
                "cell_line_name": rows.iloc[0]["cell_line_name"],
                "pairs": pairs,
                "inhibitions": inhibitions,
            }
        ]
    )


# extract combos and groupby block
combo_therapy = drug_combos.loc[
    (drug_combos["conc_c"] != 0.0) & (drug_combos["conc_r"] != 0.0)
]
combo_therapy_by_block = combo_therapy.groupby("block_id").apply(
    get_concentration_pairs
)
combo_therapy_by_block = combo_therapy_by_block.reset_index()

# connect with synergy scores
combo_therapy_by_block = combo_therapy_by_block.merge(
    combo_scores[
        [
            "block_id",
            "css_ri",
            "synergy_zip",
            "synergy_bliss",
            "synergy_loewe",
            "synergy_hsa",
            "S",
        ]
    ]
)

combo_therapy_by_block = combo_therapy_by_block[
    [
        "block_id",
        "cell_line_name",
        "drug_row_recover_id",
        "drug_row_smiles",
        "drug_col_recover_id",
        "drug_col_smiles",
        "pairs",
        "inhibitions",
        "css_ri",
        "synergy_zip",
        "synergy_bliss",
        "synergy_loewe",
        "synergy_hsa",
        "S",
    ]
]
combo_therapy_by_block = combo_therapy_by_block.rename(
    columns={"pairs": "concentration_pairs"}
)
combo_therapy_by_block.to_json(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/drug_combos/combos.json"
)

# clean summary data
summary_data = summary_data[
    ["block_id", "detection_technology", "time_h", "study_name"]
]
summary_data = summary_data.merge(
    combo_therapy_by_block[["block_id", "cell_line_name"]]
)

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

summary_data = summary_data.merge(
    mono_therapy_combined[
        ["block_id", "mono_row_measurements", "mono_col_measurements"]
    ]
)
summary_data = summary_data.merge(
    combo_therapy_by_block[["block_id", "combo_measurements"]]
)
