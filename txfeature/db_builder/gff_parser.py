"""
Takes gff annotation files as input and parses the entries and returns a dictionary.
"""


def gff_table(file):
    # Initializing the gff annotation data structure which will be used within the pipeline
    gff_ds = {'table': {}, 'tx_attr': {}}

    # Iterate through gff3 input and store information
    for gff_line in file:

        # Skip all non-data containing lines
        if gff_line[0][0] == "#":
            continue

        # The following describes the parsing instructions to create the gff data structure
        # gff is parsed and each line is stored into gff_ds under the associated gene
        # transcript to gene association is stored under tx2gene
        attr = dict(item.split('=') for item in gff_line[8].split(';'))
        gene_id = attr['gene_id']
        if gff_line[2] == 'transcript':
            gff_ds['tx_attr'][attr['transcript_id']] = {'gene_id': gene_id,
                                                        'chrom': gff_line[0],
                                                        'strand': gff_line[6],
                                                        'gene_name': attr['gene_name'],
                                                        'tx_type': attr['transcript_type']}
        if gene_id not in gff_ds['table'].keys():
            gff_ds['table'][gene_id] = []
        # Entry of gff data line into gff_ds['table']
        gff_ds['table'][gene_id].append(gff_line)

    return gff_ds
