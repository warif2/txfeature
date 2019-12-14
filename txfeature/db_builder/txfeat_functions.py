"""
contains functions for the purpose of txfeat database construction and cli functions
"""
from itertools import islice

# Function for looking up transcript features within the gff_table
def tx2gff_lookup(gff, tx):
    from .tx_classes import GffReadEntry
    entry_matches = []
    gene = gff['tx_attr'][tx]['gene_id']
    for item in gff['table'][gene]:
        gff_entry = GffReadEntry(item)
        if tx == gff_entry.transcript_id:
            entry_matches.append(item)
        else:
            continue
    return entry_matches

def chunks(data, size=10000):
    it = iter(data)
    for i in range(0, len(data), size):
        yield {i:data[i] for i in islice(it, size)}

def region_coord_aggregate(exon_dict):
    agg_coord = []
    for exon, coord in exon_dict.items():
        if coord['start'] < coord['stop']:
            agg_coord = agg_coord + list(range(coord['start'],coord['stop'] + 1))
        else:
            agg_coord = agg_coord + list(range(coord['stop'], coord['start'] + 1))
    return agg_coord