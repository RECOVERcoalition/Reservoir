import reservoir as rsv
import pandas as pd
import pickle

# get biogrid data, keep only SARS Cov-2 proteins and extract synonyms
covid_proteins_biogrid = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER
    + "/raw/proteins/covid/BIOGRID-PROJECT-covid19_coronavirus_project-GENES-4.2.191.projectindex.txt",
    sep="\t",
)
covid_proteins_biogrid = covid_proteins_biogrid.loc[
    covid_proteins_biogrid["ORGANISM NAME"]
    == "Severe acute respiratory syndrome coronavirus 2"
]
covid_proteins_biogrid["SYNONYMS"] = covid_proteins_biogrid["SYNONYMS"].apply(
    lambda synonyms: synonyms.split("|")
)
covid_proteins_biogrid = covid_proteins_biogrid.explode("SYNONYMS")
covid_proteins_biogrid["OFFICIAL SYMBOL"] = covid_proteins_biogrid[
    "OFFICIAL SYMBOL"
].apply(lambda symbol: "SARS-CoV2 " + symbol)

# create alias mapper
covid_mapper = {}
for i, row in covid_proteins_biogrid.iterrows():
    covid_mapper[row["SYNONYMS"].lower()] = row["OFFICIAL SYMBOL"]
    covid_mapper[row["OFFICIAL SYMBOL"].lower()] = row["OFFICIAL SYMBOL"]

# have to introduce some of them manually unfortunately
covid_mapper["Replicase polyprotein 1a".lower()] = "SARS-CoV2 ORF1a"
covid_mapper["Replicase polyprotein 1ab".lower()] = "SARS-CoV2 ORF1ab"
covid_mapper["Host translation inhibitor nsp1".lower()] = "SARS-CoV2 nsp1"
covid_mapper["Non-structural protein 2".lower()] = "SARS-CoV2 nsp2"
covid_mapper["Non-structural protein 3".lower()] = "SARS-CoV2 nsp3"
covid_mapper["Non-structural protein 4".lower()] = "SARS-CoV2 nsp4"
covid_mapper["Non-structural protein 6".lower()] = "SARS-CoV2 nsp6"
covid_mapper["Non-structural protein 7".lower()] = "SARS-CoV2 nsp7"
covid_mapper["Non-structural protein 8".lower()] = "SARS-CoV2 nsp8"
covid_mapper["Non-structural protein 9".lower()] = "SARS-CoV2 nsp9"
covid_mapper["Non-structural protein 10".lower()] = "SARS-CoV2 nsp10"
covid_mapper["Non-structural protein 11".lower()] = "SARS-CoV2 nsp11"
covid_mapper["2'-O-methyltransferase".lower()] = "SARS-CoV2 nsp16"
covid_mapper["Proofreading exoribonuclease".lower()] = "SARS-CoV2 nsp14"
covid_mapper["ORF10".lower()] = "SARS-CoV2 ORF10"

# load uniprot data
import json

with open(
    rsv.RESERVOIR_DATA_FOLDER
    + "/raw/proteins/covid/uniprot-download_true_format_json_query__2A_20AND_20_28other_organis-2020.10.29-22.29.57.09.json"
) as f:
    data = json.load(f)

# connect the uniprot proteins to the indentifiers in biogrid
parsed_covid_proteins = []
for protein in data["results"]:
    # get the seqence and uniprot if for the proteins
    sequence = protein["sequence"]["value"]
    uniprot_id = protein["primaryAccession"]

    # some proteins, like the replicase are split into multiple proteins. these are stored in features in the uniprot file
    if "features" in protein:
        for feature in protein["features"]:
            if feature["description"].lower() in covid_mapper:
                covid_protein = covid_mapper[feature["description"].lower()]
                start = int(feature["location"]["start"]["value"]) - 1
                end = int(feature["location"]["end"]["value"]) - 1
                parsed_covid_proteins.append(
                    {
                        "covid_protein": covid_protein,
                        "sequence": sequence[start : end + 1],
                        "uniprot_id": uniprot_id,
                    }
                )
            elif (
                "featureId" in feature and feature["featureId"].lower() in covid_mapper
            ):
                covid_protein = covid_mapper[feature["featureId"].lower()]
                start = int(feature["location"]["start"]["value"]) - 1
                end = int(feature["location"]["end"]["value"]) - 1
                parsed_covid_proteins.append(
                    {
                        "covid_protein": covid_protein,
                        "sequence": sequence[start : end + 1],
                        "uniprot_id": uniprot_id,
                    }
                )
    else:
        parsed_covid_proteins.append(
            {
                "covid_protein": covid_mapper[
                    protein["genes"][0]["geneName"]["value"].lower()
                ],
                "sequence": protein["sequence"]["value"],
                "uniprot_id": uniprot_id,
            }
        )

# clean up
parsed_covid_proteins = pd.DataFrame(parsed_covid_proteins)
parsed_covid_proteins = parsed_covid_proteins.drop_duplicates(
    "covid_protein"
).sort_values("covid_protein")

# export
parsed_covid_proteins.to_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/proteins/covid_proteins.csv", index=False
)
with open(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/proteins/covid_protein_mapper.pickle", "wb"
) as f:
    pickle.dump(covid_mapper, f)
