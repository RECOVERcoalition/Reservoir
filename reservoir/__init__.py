import reservoir

RESERVOIR_DATA_FOLDER = reservoir.__path__[0] + "/data"
from reservoir.api import *
from reservoir.parsers.drugs.molecule_mapper import map_drugs
from reservoir.parsers.proteins.protein_mapper import (
    hgnc_normalize,
    map_covid_protein,
)
from reservoir.parsers.drugs.molecule_tools import *
