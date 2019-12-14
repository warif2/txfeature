"""
transcript_assembly.py assembles given transcripts using the provided annotation and fasta file.
"""

import csv
import copy
import os
from pybedtools import BedTool
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from .tx_classes import GffReadEntry


def build(transcript, annot, txi_dict, fasta_path, tmp_dir=''):
    # Create temp_dir
    if tmp_dir == '':
        tmp_dir = './_tmp_txfeature'
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # Initializing return variables
    transcript_status = {'five_prime_UTR': '', 'stop_codon_redefined_as_selenocysteine': '', 'exon': '',
                         'stop_codon': '', 'CDS': '', 'three_prime_UTR': '', 'start_codon': ''}

    # tx associated variables needed for assembly
    chrom = txi_dict['chrom']
    strnd = txi_dict['strand']
    gene_id = txi_dict['gene_id']
    tx_type = txi_dict['tx_type']

    # Find all features of transcript from gff_table and store it into tran_features dictionary
    tx_annot = {'five_prime_UTR': {}, 'stop_codon_redefined_as_selenocysteine': {}, 'exon': {},
                     'stop_codon': {}, 'CDS': {}, 'three_prime_UTR': {}, 'start_codon': {}}
    seleno_index = 0
    for item in annot:
        gff_feature = GffReadEntry(item)
        ftype = gff_feature.entry_type
        if ftype in ['five_prime_UTR', 'exon', 'CDS', 'three_prime_UTR', 'stop_codon', 'start_codon']:
            tx_annot[ftype][gff_feature.exon_number] = {'start': gff_feature.start_coord,
                                                             'stop': gff_feature.stop_coord}
        if ftype in ['stop_codon_redefined_as_selenocysteine']:
            tx_annot[ftype][seleno_index] = {'start': gff_feature.start_coord, 'stop': gff_feature.stop_coord}
            seleno_index -= 1

    # Complete transcript_status
    for ftype in transcript_status.keys():
        if len(tx_annot[ftype]) == 0:
            transcript_status[ftype] = 'not_defined'
        else:
            transcript_status[ftype] = 'defined'

    # Assemble transcript using tran_features
    # Create BED file of from the 'exon' features in tran_features
    # name = ex_num:ex_start_coord:ex_stop_coord, this will be used when BED file is passed into pybedtools_getfasta
    with open(tmp_dir + '/' + transcript + '_temp.bed', 'w') as temp_bed:
        bed_out = csv.writer(temp_bed, delimiter='\t')
        for exon_number in range(1, max(tx_annot['exon'].keys()) + 1):
            exon_coord = tx_annot['exon'][exon_number]
            name = '{ex_num}:{start}:{stop}'.format(ex_num=exon_number, start=exon_coord['start'],
                                                    stop=exon_coord['stop'])
            bed_out.writerow([chrom, exon_coord['start'] - 1, exon_coord['stop'], name, 0, strnd])
    temp_bed.close()

    # Pass BED file into pybedtools getfasta and extract transcript sequence
    transcript_bed = BedTool(tmp_dir + '/' + transcript + '_temp.bed')
    fasta_file = BedTool(fasta_path)
    transcript_fasta = transcript_bed.sequence(fi=fasta_file, name=True, s=True,
                                               fo=tmp_dir + '/' + transcript + '_temp.fa')

    # Parse transcript fasta and create g2iloc table
    g2iloc = []
    full_seq = ''
    with open(transcript_fasta.seqfn) as transcript_seq:
        for fasta_entry in transcript_seq:
            if fasta_entry[0] == '>':
                seq_info = fasta_entry.strip('>').strip('\n').replace('(', ':').split(':')
                ex_seq = next(transcript_seq).strip('\n')
                full_seq += ex_seq
                if strnd == '+':
                    seq_start = int(seq_info[1])
                if strnd == '-':
                    seq_start = int(seq_info[2])
                for base in range(len(ex_seq)):
                    g2iloc.append(seq_start)
                    if strnd == '+':
                        seq_start += 1
                    if strnd == '-':
                        seq_start -= 1
    g2iloc = tuple(g2iloc)
    transcript_seq.close()

    # Setup output dictionary
    output = {'tx_id': transcript,
              'gene_id': gene_id,
              'tx_type': tx_type,
              'chrom': chrom,
              'strand': strnd,
              'tx_status': transcript_status,
              'tx_seq': full_seq,
              'num_of_exons': max(tx_annot['exon'].keys()),
              'g2iloc': g2iloc,
              'tx_annot': tx_annot}
    return output


def translate_seq(transcript_seq):
    # Check if CDS, start_codon and stop_codon are all defined
    transcript_status = TranscriptRead.Status(transcript_seq['transcript_status'])
    if (transcript_status.status_cds, transcript_status.status_start, transcript_status.status_stop) == \
            ('defined', 'defined', 'defined'):
        pass
    else:
        return ['CDS not well defined!']

    # Load 'n_seq' into transcript_ds using exfeat_classes.transcript_read. Use mrna_region method to obtain CDS
    # sequence
    transcript_ds = TranscriptRead.Sequence(transcript_seq['n_seq'])
    cds_index = transcript_ds.index(region_type='CDS')
    cds_seq = transcript_ds.mrna_region('CDS')

    # Translate cds_seq
    coding_dna = Seq(cds_seq, generic_dna)
    protein_seq = coding_dna.translate()

    # Create data structure for amino acid sequence of translated protein and store it in protein_ds
    protein_ds = []

    # Iterate through each amino acid of protein_seq
    for index, amino_acid in enumerate(protein_seq, 1):
        coords = []
        exons = []

        # Iterate through each base of the associated codon and extract exon and coordinate information
        for transcript_index in range(cds_index[0] + ((index - 1) * 3), cds_index[0] + (index * 3)):
            cds_base = TranscriptRead.Sequence.Base(transcript_ds.transcript_seq[transcript_index])
            coords.append(cds_base.base_coord)
            if cds_base.exon_number not in exons:
                exons.append(cds_base.exon_number)

        # Setup return array
        protein_ds.append([index,  # index of the amino acid in protein sequence
                           amino_acid,  # amino acid letter
                           cds_seq[(index - 1) * 3:index * 3],  # codon
                           copy.deepcopy(coords),  # coord of transcript coding amino acid
                           copy.deepcopy(exons)])  # exons of transcript

    return protein_ds


if __name__ == '__main__':
    pass
