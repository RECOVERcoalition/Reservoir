import reservoir as rsv
import pandas as pd

"""
Gordon Krogan
"""
krogan_data = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/raw/ppi/covid_gordon_krogan_ppi.csv"
)

# normalize protein ids
krogan_data["human_gene_hgnc_id"] = krogan_data["PreyGene"].apply(rsv.hgnc_normalize)
krogan_data = krogan_data[["Bait", "human_gene_hgnc_id"]].rename(
    columns={"Bait": "covid_protein"}
)

# fix some of the names
krogan_data.loc[
    krogan_data["covid_protein"] == "SARS-CoV2 Spike", "covid_protein"
] = "SARS-CoV2 S"
krogan_data.loc[
    krogan_data["covid_protein"] == "SARS-CoV2 nsp5_C145A", "covid_protein"
] = "SARS-CoV2 nsp5"

assert krogan_data.loc[krogan_data["human_gene_hgnc_id"].isna()].shape[0] == 0

# map to covid protein id
krogan_data["covid_protein"] = krogan_data["covid_protein"].apply(rsv.map_covid_protein)
krogan_data.to_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/ppi/covid_krogan_ppi.csv", index=False
)


"""
Li Liang
"""
liang_data = pd.read_csv(rsv.RESERVOIR_DATA_FOLDER + "/raw/ppi/covid_li_liang_ppi.csv")

# normalize human names
liang_data["human_gene_hgnc_id"] = liang_data["Gene name"].apply(rsv.hgnc_normalize)
liang_data = liang_data.loc[~liang_data["human_gene_hgnc_id"].isna()]

# fix covid proteins
liang_data = liang_data[["Bait", "human_gene_hgnc_id"]].rename(
    columns={"Bait": "covid_protein"}
)
liang_data["covid_protein"] = liang_data["covid_protein"].apply(
    lambda name: "SARS-CoV2 " + name
)
liang_data.loc[
    liang_data["covid_protein"] == "SARS-CoV2 nsp3-c", "covid_protein"
] = "SARS-CoV2 nsp3"
liang_data.loc[
    liang_data["covid_protein"] == "SARS-CoV2 nsp3-n", "covid_protein"
] = "SARS-CoV2 nsp3"

# map to covid protein id
liang_data["covid_protein"] = liang_data["covid_protein"].apply(rsv.map_covid_protein)

liang_data.to_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/ppi/covid_liang_ppi.csv", index=False
)
