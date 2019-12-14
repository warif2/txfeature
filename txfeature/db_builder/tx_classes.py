"""
exfeat_classes.py contain the major classes used within the ExFeat Pipeline.
"""
import txfeature.db_builder.txfeat_functions as txfeat_func


class GffReadEntry:

    def __init__(self, entry):

        # General feature information
        self.chr = entry[0]
        self.source = entry[1]
        self.entry_type = entry[2]
        self.start_coord = int(entry[3])
        self.stop_coord = int(entry[4])
        self.strand = entry[6]

        # Converting attribute column into dictionary
        entry_attr = {}
        for tag in entry[8].split(';'):
            tag_id = tag.split('=')[0]
            tag_value = tag.split('=')[1]
            entry_attr[tag_id] = tag_value

        # General Gene Attributes
        self.gene_type = entry_attr['gene_type']
        self.gene_name = entry_attr['gene_name']
        self.gene_id = entry_attr['gene_id']
        if 'gene_status' in entry_attr.keys():
            self.gene_status = entry_attr['gene_status']
        else:
            self.gene_status = "NA"
        if 'havana_gene' in entry_attr.keys():
            self.havana_gene = entry_attr['havana_gene']
        else:
            self.havana_gene = "NA"
        self.level = entry_attr['level']

        # Transcript, Exon, CDS, 5UTR, 3UTR, Start_codon, Stop_codon, and Seleno_stop_codon Attributes
        self.transcript_id = "NA"
        if self.entry_type in ['transcript', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'start_codon',
                               'stop_codon', 'stop_codon_redefined_as_selenocysteine']:
            self.transcript_id = entry_attr['transcript_id']
            self.parent = entry_attr['Parent']
            self.transcript_type = entry_attr['transcript_type']
            if self.transcript_type == 'protein_coding':
                self.protein_id = entry_attr['protein_id']
            else:
                self.protein_id = "NA"
            self.transcript_name = entry_attr['transcript_name']
            if 'transcript_status' in entry_attr.keys():
                self.transcript_status = entry_attr['transcript_status']
            else:
                self.transcript_status = "NA"
            if 'havana_transcript' in entry_attr.keys():
                self.havana_transcript = entry_attr['havana_transcript']
            else:
                self.havana_transcript = 'NA'
            if 'transcript_support_level' in entry_attr.keys():
                self.transupp_level = entry_attr['transcript_support_level']
            else:
                self.transupp_level = "NA"
            if 'tag' in entry_attr.keys():
                self.tag = entry_attr['tag']
            else:
                self.tag = "NA"

            if self.entry_type in ['exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'start_codon', 'stop_codon',
                                   'stop_codon_redefined_as_selenocysteine']:
                if self.entry_type not in ['stop_codon_redefined_as_selenocysteine']:
                    self.exon_number = int(entry_attr['exon_number'])
                    self.exon_id = entry_attr['exon_id']
                if self.entry_type in ['CDS', 'five_prime_UTR', 'three_prime_UTR', 'start_codon', 'stop_codon',
                                       'stop_codon_redefined_as_selenocysteine']:
                    if entry[7] == '.':
                        self.phase = entry[7]
                    else:
                        self.phase = int(entry[7])


class TxRead:

    def __init__(self, transcript):

        # Transcript associated variables
        self.tx_id = transcript['tx_id']
        self.gene_id = transcript['gene_id']
        self.tx_type = transcript['tx_type']
        self.sequence = transcript['tx_seq']
        self.num_exons = transcript['num_of_exons']
        self.g2iloc = transcript['g2iloc']
        self.tx_status = transcript['tx_status']
        self.tx_annot = transcript['tx_annot']
        self.chrom = transcript['chrom']
        self.strand = transcript['strand']
        if self.tx_status['start_codon'] == 'defined':
            self.start_codon_exon = list(self.tx_annot['start_codon'].keys())[0]
            self.start_codon_coord = self.chrom + ':' + \
                                     str(self.tx_annot['start_codon'][self.start_codon_exon]['start']) + '-' + \
                                     str(self.tx_annot['start_codon'][self.start_codon_exon]['stop'])
        else:
            self.start_codon_exon = 'NA'
            self.start_codon_coord = 'NA'
        if self.tx_status['stop_codon'] == 'defined':
            self.stop_codon_exon = list(self.tx_annot['stop_codon'].keys())[0]
            self.stop_codon_coord = self.chrom + ':' + \
                                    str(self.tx_annot['stop_codon'][self.stop_codon_exon]['start']) + '-' + \
                                    str(self.tx_annot['stop_codon'][self.stop_codon_exon]['stop'])
        else:
            self.stop_codon_exon = 'NA'
            self.stop_codon_coord = 'NA'

    # Method to retreive region information from genomic_coord
    def coord_to_region(self, gen_coord):
        # setting up return variable
        region = []
        # read in genomic coordinate
        qchrom = gen_coord.split(':')[0]
        qcoord = gen_coord.split(':')[1].split('-')
        if len(qcoord) == 1:
            qcoord = {int(qcoord[0])}
        elif len(qcoord) == 2:
            qcoord = set(list(range(int(qcoord[0]), int(qcoord[1]) + 1)))
        else:
            return region  # invalid coordinate, return empty
        # check overlap
        if qchrom == self.chrom:
            for region_type in self.tx_annot.keys():
                region_coord = set(txfeat_func.region_coord_aggregate(self.tx_annot[region_type]))
                if len(qcoord.intersection(region_coord)) > 0:
                    region.append(region_type)
        return region

    # Method to extract sequence of specified mRNA region, ex 'five_prime_UTR'
    def get_sequence(self, region_type, query):
        # setting return variable
        sequence = ''
        # type for mRNA regions
        if region_type == 'mrna_region':
            if query in self.tx_annot.keys():
                region_coord = txfeat_func.region_coord_aggregate(self.tx_annot[query])
            else:
                return sequence
            if self.strand == '+':
                start = min(region_coord)
                end = max(region_coord)
            else:
                start = max(region_coord)
                end = min(region_coord)
            seqi_start = self.g2iloc.index(start)
            seqi_end = self.g2iloc.index(end)

        elif region_type == 'exon':
            if query in self.tx_annot['exon'].keys():
                ex_coord = self.tx_annot['exon'][query]
                if self.strand == '+':
                    seqi_start = self.g2iloc.index(ex_coord['start'])
                    seqi_end = self.g2iloc.index(ex_coord['stop'])
                else:
                    seqi_start = self.g2iloc.index(ex_coord['stop'])
                    seqi_end = self.g2iloc.index(ex_coord['start'])
            else:
                return sequence

        elif region_type == 'coordinates':
            qchrom = query.split(':')[0]
            if qchrom != self.chrom:
                return ''
            qcoord = query.split(':')[1].split('-')
            if len(qcoord) == 1:
                qcoord = int(qcoord[0])
                if qcoord in self.g2iloc:
                    index = self.g2iloc.index(qcoord)
                    return self.sequence[index]
                else:
                    return ''  # coordinate not in transcript
            elif len(qcoord) == 2:
                start = int(qcoord[0])
                end = int(qcoord[1])
                if start in self.g2iloc and end in self.g2iloc:
                    if self.strand == '+':
                        seqi_start = self.g2iloc.index(start)
                        seqi_end = self.g2iloc.index(end)
                    else:
                        seqi_start = self.g2iloc.index(end)
                        seqi_end = self.g2iloc.index(start)
                else:
                    return ''  # coordinates not in transcript
            else:
                return ''  # invalid coordinate, return empty

        elif region_type == 'coordinates_overlap':
            qchrom = query.split(':')[0]
            if qchrom != self.chrom:
                return ''
            qcoord = query.split(':')[1].split('-')
            if len(qcoord) == 2:
                start = int(qcoord[0])
                end = int(qcoord[1])
                qspan = set(list(range(start, end + 1)))
                ovlp = qspan.intersection(set(self.g2iloc))
                if self.strand == '+':
                    seqi_start = self.g2iloc.index(min(ovlp))
                    seqi_end = self.g2iloc.index(max(ovlp))
                else:
                    seqi_start = self.g2iloc.index(max(ovlp))
                    seqi_end = self.g2iloc.index(min(ovlp))
            else:
                return ''  # invalid coordinate, return empty

        else:
            return ''
        return self.sequence[seqi_start:(seqi_end + 1)]

    def length(self, region_type, query):
        # type for mRNA regions
        if region_type == 'mrna_region':
            if query in self.tx_annot.keys():
                region_coord = txfeat_func.region_coord_aggregate(self.tx_annot[query])
                return len(region_coord)
            else:
                return 0

        elif region_type == 'exon':
            if query in self.tx_annot['exon'].keys():
                ex_coord = self.tx_annot['exon'][query]
                return abs(ex_coord['start'] - ex_coord['stop']) + 1
            else:
                return 0


class ProteinRead:
    # Protein sequence
    class Sequence:
        def __init__(self, protein_sequence):
            self.protein_seq = protein_sequence

        # Method to extract sequence of specified mRNA region, ex 'five_prime_UTR'
        def exon_number(self, exons):
            ret_seq = ''
            for seq in self.protein_seq:
                amino = self.AminoAcid(seq)
                if len(set(exons).intersection(amino.exon)) > 0:
                    ret_seq += amino.aa_letter
            if len(ret_seq) == 0:
                return 'not_defined'
            return ret_seq

        # Protein Amino Acid
        class AminoAcid:
            def __init__(self, protein_aa):
                self.index = protein_aa[0]
                self.aa_letter = protein_aa[1]
                self.codon = protein_aa[2]
                self.transcript_coord = protein_aa[3]
                self.exon = protein_aa[4]
