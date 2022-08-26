import reservoir as rsv
import pandas as pd
import networkx as nx
import itertools
import pickle


def clean_cell_line_name(cell_line_name):
    """
    Create versions that don't have digits and numbers
    A-375 and A 375 => A375
    """
    clean_name = str(cell_line_name).lower()

    if "ncih" in clean_name and "_" in clean_name:
        clean_name = clean_name.split("_")[0]

    for ch in [
        " ",
        "(",
        ",",
        "\\",
        "{",
        ".",
        "~",
        "?",
        "]",
        "*",
        "^",
        ")",
        " ",
        "[",
        "#",
        "/",
        "_",
        "-",
        "'",
        ":",
        ";",
        "+",
        "}",
    ]:
        clean_name = clean_name.replace(ch, "")

    return clean_name


if __name__ == "__main__":
    cd = pd.read_csv(
        rsv.RESERVOIR_DATA_FOLDER + "/raw/cell_lines/cell_line_dictionary_synonyms.tab",
        sep="\t",
    )
    cd["symbol"] = cd["symbol"].str.lower()

    # connect cell lines based on clean name
    cd["clean_symbol"] = cd["symbol"].apply(clean_cell_line_name)
    g = nx.Graph()
    g.add_nodes_from(cd["symbol"])

    grouped_symbols = cd.groupby("clean_symbol").agg({"symbol": set})
    grouped_symbols = grouped_symbols.loc[grouped_symbols["symbol"].apply(len) > 1]
    for symbols in grouped_symbols["symbol"]:
        for s1, s2 in itertools.combinations(symbols, 2):
            g.add_edge(s1, s2)

    # compress the cell lines based on connected components
    cell_lines_compressed = []
    for i, group in enumerate(nx.connected_components(g)):
        for cell_line in group:
            cell_lines_compressed.append(
                {
                    "cell_line": cell_line,
                    "cell_line_clean": clean_cell_line_name(cell_line),
                    "cell_line_id": str(i),
                }
            )

    cell_lines_compressed = pd.DataFrame(cell_lines_compressed)
    cell_lines_compressed = cell_lines_compressed[
        ["cell_line_clean", "cell_line_id"]
    ].drop_duplicates()
    cd = cell_lines_compressed.set_index("cell_line_clean").to_dict()["cell_line_id"]

    with open(rsv.RESERVOIR_DATA_FOLDER + "/parsed/cell_lines/cd.pickle", "wb") as f:
        pickle.dump(cd, f)
