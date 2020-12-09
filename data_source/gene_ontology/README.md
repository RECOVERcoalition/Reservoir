# Gene ontology
The Gene Ontology (GO) knowledge base is the world’s largest source of information on the functions of genes
The functions are annotated to genes using GO terms in the format GO:1904469. The ontology is hierarchical, and annotations usually happen at leaf nodes so they need to get propagated

## Structure
The ontology has the following raw files:
* `goa_human.gaf`: file with the annotations Gene->GO Term
* `go.obo`: file with the hierarchical ontology, basically an directed acyclical graph
