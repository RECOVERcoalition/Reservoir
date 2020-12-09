import pandas as pd
import reservoir as rsv
from reservoir.parsers.proteins import protein_mapper
import ppi_helper

# load data
biogrid_mv = pd.read_csv(
    rsv.RESERVOIR + "/raw/ppi/BIOGRID-MV-Physical-4.0.189.tab3.txt", sep="\t"
)
biogrid = pd.read_csv(
    rsv.RESERVOIR + "/raw/ppi/BIOGRID-ORGANISM-Homo_sapiens-4.0.189.tab3.txt", sep="\t",
)

# keep only human - human interactions
biogrid = biogrid.loc[
    (biogrid["Organism ID Interactor A"] == 9606)
    & (biogrid["Organism ID Interactor B"] == 9606),
]
biogrid_mv = biogrid_mv.loc[
    (biogrid_mv["Organism ID Interactor A"] == 9606)
    & (biogrid_mv["Organism ID Interactor B"] == 9606)
]

# normalize protein ids and keep only valid ones
biogrid["gene_1_hgnc_id"] = biogrid["SWISS-PROT Accessions Interactor A"].apply(
    protein_mapper.uniprot_to_hgnc
)
biogrid["gene_2_hgnc_id"] = biogrid["SWISS-PROT Accessions Interactor B"].apply(
    protein_mapper.uniprot_to_hgnc
)
biogrid = biogrid.loc[
    ~pd.isnull(biogrid["gene_1_hgnc_id"]) & ~pd.isnull(biogrid["gene_2_hgnc_id"])
]

biogrid_mv["gene_1_hgnc_id"] = biogrid["SWISS-PROT Accessions Interactor A"].apply(
    protein_mapper.uniprot_to_hgnc
)
biogrid_mv["gene_2_hgnc_id"] = biogrid["SWISS-PROT Accessions Interactor B"].apply(
    protein_mapper.uniprot_to_hgnc
)
biogrid_mv = biogrid_mv.loc[
    ~pd.isnull(biogrid_mv["gene_1_hgnc_id"]) & ~pd.isnull(biogrid_mv["gene_2_hgnc_id"])
]

# retain only relevant coumns
biogrid = biogrid[
    [
        "gene_1_hgnc_id",
        "gene_2_hgnc_id",
        "Experimental System",
        "Experimental System Type",
        "Author",
        "Publication Source",
    ]
]
biogrid_mv = biogrid_mv[["gene_1_hgnc_id", "gene_2_hgnc_id"]]

# remove duplicates for the mv file only
biogrid_mv = ppi_helper.remove_duplicate_edges(biogrid_mv)

biogrid.to_csv(rsv.RESERVOIR + "/parsed/ppi/biogrid.csv", index=False)
biogrid_mv.to_csv(rsv.RESERVOIR + "/parsed/ppi/biogrid_mv.csv", index=False)
