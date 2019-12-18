"""
contains functions to determine properties associated to transcript sequences
"""

import re
from txfeature.db_builder import utils


def au_element(seq):
    """
    Takes input string and calculates 1) 'AUUUA' ppentamer count 2) number of AU elements 3) fraction of sequence
    consisting of AU element and 4) longest AU element in sequence. AU element is considered a sequence of A or U longer
    than 5 bases.
    :param seq: input mRNA sequence as string
    :returns: dict of the calculated values. dict keys: 'au_pentamer', 'au_num', 'au_fraction', 'au_longest'
    """
    # convert T to U
    seq = seq.replace('T', 'U')

    # pentamer count
    au_pentamer_count = seq.count('AUUUA')

    # search au elements
    aurich_re = re.compile(r'[AU]{5,}')
    au_search = aurich_re.findall(seq)

    # calculate values
    num_au_elements = len(au_search)
    au_fraction = len("".join(au_search)) / float(len(seq))
    if num_au_elements > 0:
        longest_au_element = max(map(len, au_search))
    else:
        longest_au_element = 0

    return {'au_pentamer': au_pentamer_count,
            'au_num': num_au_elements,
            'au_fraction': au_fraction,
            'au_longest': longest_au_element}


def gc_content(seq):
    """
    Calculates % GC content of sequence ignoring any N present in sequence.
    :param seq: input mRNA sequence as string
    :return: percent gc as float
    """
    # count G and C in seq
    num_g = seq.count('G')
    num_c = seq.count('C')
    num_n = seq.count('N')

    # calculate gc content
    percent_gc = (num_g + num_c) / float(len(seq) - num_n) if len(seq) - num_n > 0 else 0

    return percent_gc


def rnafold_energy(sequence, command):
    if len(sequence) < 1:
        return {'mfe': 0, 'ensemble': 0, 'centroid': 0, 'mea': 0}

    out_array = []
    output = utils.stdout_from_command("echo %s | %s" % (sequence, command))
    for lines in output:
        out_array.append(str(lines, 'utf-8'))

    # second line has the MFE
    mfe = float(re.sub('[()]', '', out_array[1].split()[-1]))

    # third line has ensemble energy
    ensemble = out_array[2].split()[-1].strip('[]')
    #ensemble = float(re.sub(r'([*?])', '', out_array[2].split()[-1]))

    # fourth line has the centroid energy
    centroid = float(re.sub('[{}]', '', out_array[3].split()[-2]))

    # fifth line has the MEA
    mea = float(re.sub('[{}]', '', out_array[4].split()[-2]))

    return {'mfe': mfe, 'ensemble': ensemble, 'centroid': centroid, 'mea': mea}


def rnalfold_energy(sequence, command):
    if len(sequence) < 1:
        return 0

    out_array = []
    output = utils.stdout_from_command("echo %s | %s" % (sequence, command))
    for lines in output:
        out_array.append(str(lines, 'utf-8'))
    pattern = re.compile(r'(?<=\().*?(?=\))')
    energy = float(pattern.findall(out_array[-1])[0])
    return energy


def kozac_score(txread):
    if txread.tx_status['five_prime_UTR'] == 'defined' and txread.tx_status['start_codon'] == 'defined':
        if txread.length('mrna_region', 'five_prime_UTR') > 6 and txread.length('mrna_region', 'CDS') > 6:
            if txread.strand == '+':
                start_zero = txread.g2iloc.index(int(txread.start_codon_coord.split(':')[1].split('-')[1]))
            else:
                start_zero = txread.g2iloc.index(int(txread.start_codon_coord.split(':')[1].split('-')[0]))
            score_1 = 0
            score_3 = 0
            pos_1 = txread.sequence[start_zero + 1]  # +3 if G
            neg_3 = txread.sequence[start_zero - 3]  # +1 if C
            neg_4 = txread.sequence[start_zero - 4]  # +1 if C
            neg_5 = txread.sequence[start_zero - 5]  # +3 if A or G
            neg_6 = txread.sequence[start_zero - 6]  # +1 if C
            neg_7 = txread.sequence[start_zero - 7]  # +1 if C
            neg_8 = txread.sequence[start_zero - 8]  # +3 if G
            for position in [neg_3, neg_4, neg_6, neg_7]:
                if position == 'C':
                    score_1 += 1
            for position in [pos_1, neg_8]:
                if position == 'G':
                    score_3 += 1
            if neg_5 == 'A' or neg_5 == 'G':
                score_3 += 1
            kozac_score = score_1 + 3*(score_3)
        else:
            return -1
    else:
        return -1

    return kozac_score


def top_score(sequence):
    # TODO Complete TOP score algo
    pass

if __name__ == '__main__':
    print(au_element('ATTUUA'))
