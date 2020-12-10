import recover_data_lake

RECOVER_DATA_FOLDER = recover_data_lake.__path__[0] + "/data"
from recover_data_lake.api import *
from recover_data_lake.parsers.drugs.molecule_mapper import map_drugs
from recover_data_lake.parsers.proteins.protein_mapper import (
    hgnc_normalize,
    map_covid_protein,
)
from recover_data_lake.parsers.drugs.molecule_tools import *
