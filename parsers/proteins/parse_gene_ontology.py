import pandas as pd
import reservoir as rsv
from reservoir.parsers.proteins import protein_mapper
import obonet
import networkx as nx

# normalize the hgnc symbol in the go annotations
go_annotations = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/raw/gene_ontology/goa_human.gaf",
    comment="!",
    sep="\t",
    names=list(map(str, range(17))),
)
go_annotations["2"] = go_annotations["2"].apply(protein_mapper.hgnc_normalize)

# keep only valid proteins
human_proteins = set(
    pd.read_csv(rsv.RESERVOIR_DATA_FOLDER + "/parsed/proteins/human_proteins.csv")[
        "gene_hgnc_id"
    ]
)
go_annotations = go_annotations.loc[
    go_annotations["2"].apply(lambda gene: gene in human_proteins)
]
go_annotations = go_annotations[["2", "4", "6"]].rename(
    columns={"2": "gene_hgnc_id", "4": "go_id", "6": "evidence"}
)

# create code mapper
codes = {}
codes["experimental evidence"] = [
    "EXP",
    "IDA",
    "IPI",
    "IMP",
    "IGI",
    "IEP",
    "HTP",
    "HDA",
    "HMP",
    "HGI",
    "HEP",
]
codes["phylogenetically inferred"] = ["IBA", "IBD", "IKR", "IRD"]
codes["computational analysis"] = ["ISS", "ISO", "ISA", "ISM", "IGC", "RCA"]
codes["author statement"] = ["TAS", "NAS"]
codes["curator statement"] = ["IC", "ND"]
codes["electronic annotation"] = ["IEA"]
code_mapper = {code: category for category in codes for code in codes[category]}
go_annotations["evidence_category"] = go_annotations["evidence"].apply(
    lambda evidence: code_mapper[evidence]
)

# create ontology graph
go_graph = obonet.read_obo(rsv.RESERVOIR_DATA_FOLDER + "/raw/gene_ontology/go.obo")

# keep only is_a and part_of edges
for (a, b, t) in list(go_graph.edges):
    if t not in ["is_a", "part_of"]:
        go_graph.remove_edge(a, b, t)

# calculate the distance from the three different roots
bp_distance = pd.DataFrame(
    [
        {
            "go_id": source,
            "distance_from_root": distance,
            "aspect": "biological_process",
        }
        for source, distance in nx.single_target_shortest_path_length(
            go_graph, target="GO:0008150"
        )
        if go_graph.nodes[source]["namespace"] == "biological_process"
    ]
)
mf_distance = pd.DataFrame(
    [
        {
            "go_id": source,
            "distance_from_root": distance,
            "aspect": "molecular_function",
        }
        for source, distance in nx.single_target_shortest_path_length(
            go_graph, target="GO:0003674"
        )
        if go_graph.nodes[source]["namespace"] == "molecular_function"
    ]
)
cc_distance = pd.DataFrame(
    [
        {
            "go_id": source,
            "distance_from_root": distance,
            "aspect": "cellular_component",
        }
        for source, distance in nx.single_target_shortest_path_length(
            go_graph, target="GO:0005575"
        )
        if go_graph.nodes[source]["namespace"] == "cellular_component"
    ]
)
assert len(set(bp_distance["go_id"]).intersection(set(mf_distance["go_id"]))) == 0
assert len(set(mf_distance["go_id"]).intersection(set(cc_distance["go_id"]))) == 0
assert len(set(bp_distance["go_id"]).intersection(set(cc_distance["go_id"]))) == 0

distance_df = pd.concat([bp_distance, mf_distance, cc_distance])

# add descendant annotations
descendant_annotations = []
for i, leaf_row in go_annotations.iterrows():
    if leaf_row["go_id"] in go_graph.nodes:
        for parent_go_id in nx.descendants(go_graph, leaf_row["go_id"]):
            descendant_annotations.append(
                {
                    "gene_hgnc_id": leaf_row["gene_hgnc_id"],
                    "go_id": parent_go_id,
                    "evidence": leaf_row["evidence"],
                    "evidence_category": leaf_row["evidence_category"],
                }
            )
# combine with leaf annotations
descendant_annotations = pd.DataFrame(descendant_annotations)
full_annotations = pd.concat([go_annotations, descendant_annotations]).drop_duplicates()

# add distances and aspects
full_annotations = full_annotations.merge(distance_df, how="left")
full_annotations["distance_from_root"] = full_annotations["distance_from_root"].fillna(
    -1
)
full_annotations = full_annotations.astype({"distance_from_root": "int32"})

# export
full_annotations.to_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/proteins/gene_ontology.csv", index=False
)
