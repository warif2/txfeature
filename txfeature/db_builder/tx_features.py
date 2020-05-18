"""
Takes assembled transcripts and returns dictionary of transcript features.
"""

import logging
import time

from txfeature.db_builder import utils, tx_classes
from txfeature.db_builder import txseq_properties as tp
from txfeature.db_builder import build_config


# def add_features(tx_assembled, tx_dict, job):
def add_features(tx_assembled, job, chunk_len, threads):
    """
    :param tx_assembled: parsed gff file using gff_parser
    :param tx_dict: list of transcript id for assembly
    :param job: specifies job number for multiprocessing
    return dict of transcript with features
    """

    # setup logger and time
    logger = logging.getLogger(__name__ + '.add_features')
    # logger.debug('Feature aggregation job %i initiated...' % job)
    initial_time = time.time()

    # setup return structure and tx list
    tx_feature = []

    # progress bar variables
    # progress = 0
    # total_tx = len(tx_dict.keys())
    stepper = 0

    # set up commands
    rnafold_command = build_config.cfg['viennarna_dir'] + build_config.cfg['rnafold_command']
    rnalfold_command = build_config.cfg['viennarna_dir'] + build_config.cfg['rnalfold_command']

    # iterate through tx_list and construct table
    for tx in tx_assembled.keys():
        # display progress of txfeat construction
        if job % threads == 0:
            utils.waiting_bar(stepper, message='Chunk %i / %i completed' % (job, chunk_len))

        # load transcript into TxRead class
        tx_read = tx_classes.TxRead(tx_assembled[tx])

        # initialize tx_feature dict
        txfeat_dict = {'tx_id': tx_read.tx_id,
                       'gene_id': tx_read.gene_id,
                       'tx_type': tx_read.tx_type}

        try:
            # get all transcript features
            txfeat_dict.update({'tx.length': len(tx_read.sequence),
                                'tx.exon_count': tx_read.num_exons,
                                'tx.gc': tp.gc_content(tx_read.sequence)})

            if tx_read.tx_type == 'protein_coding':
                txfeat_dict.update({'tx.kozac_score': tp.kozac_score(tx_read)})
            else:
                txfeat_dict.update({'tx.kozac_score': 'NA'})
        except:
            logger.debug('Job %s Error: Transcript %s features could not be determined' % (job, tx_read.tx_id))

        # get utr5 features
        try:
            if tx_read.tx_status['five_prime_UTR'] == 'defined':
                utr5_sequence = tx_read.get_sequence('mrna_region', 'five_prime_UTR')
                scan_energy = tp.rnalfold_energy(utr5_sequence, rnalfold_command)
                if len(utr5_sequence) > 50:
                    cap_energy = tp.rnafold_energy(utr5_sequence[0:50], rnafold_command)
                else:
                    cap_energy = tp.rnafold_energy(utr5_sequence, rnafold_command)
                txfeat_dict.update({'utr5.length': tx_read.length('mrna_region', 'five_prime_UTR'),
                                    'utr5.gc': tp.gc_content(utr5_sequence),
                                    'utr5.cap_structure_mfe': cap_energy['mfe'],
                                    'utr5.structure_min_scan': scan_energy})
                if len(utr5_sequence) < 1500:
                    region_energy = tp.rnafold_energy(utr5_sequence, rnafold_command)
                    txfeat_dict.update({'utr5.structure_mfe': region_energy['mfe'],
                                        'utr5.structure_mea': region_energy['mea'],
                                        'utr5.structure_centroid': region_energy['centroid']})
        except:
            logger.debug('Job %s Error: Transcript %s five_prime_UTR feature undetermined' % (job, tx_read.tx_id))

        # get cds features
        try:
            if tx_read.tx_status['CDS'] == 'defined':
                cds_sequence = tx_read.get_sequence('mrna_region', 'CDS')
                scan_energy = tp.rnalfold_energy(cds_sequence, rnalfold_command)
                txfeat_dict.update({'cds.length': tx_read.length('mrna_region', 'CDS'),
                                    'cds.gc': tp.gc_content(cds_sequence),
                                    'cds.structure_min_scan': scan_energy})
                if len(cds_sequence) < 1500:
                    region_energy = tp.rnafold_energy(cds_sequence, rnafold_command)
                    txfeat_dict.update({'cds.structure_mfe': region_energy['mfe'],
                                        'cds.structure_mea': region_energy['mea'],
                                        'cds.structure_centroid': region_energy['centroid']})
        except:
            logger.debug('Job %s Error: Transcript %s cds feature undetermined' % (job, tx_read.tx_id))

        # get utr3 features
        try:
            if tx_read.tx_status['three_prime_UTR'] == 'defined':
                utr3_sequence = tx_read.get_sequence('mrna_region', 'three_prime_UTR')
                scan_energy = tp.rnalfold_energy(utr3_sequence, rnalfold_command)
                au_element = tp.au_element(utr3_sequence)
                txfeat_dict.update({'utr3.length': tx_read.length('mrna_region', 'CDS'),
                                    'utr3.gc': tp.gc_content(utr3_sequence),
                                    'utr3.au_pentamer': au_element['au_pentamer'],
                                    'utr3.au_count': au_element['au_num'],
                                    'utr3.au_fraction': au_element['au_fraction'],
                                    'utr3.au_longest': au_element['au_longest'],
                                    'utr3.strucutre_min_scan': scan_energy})
                if len(utr3_sequence) < 1500:
                    region_energy = tp.rnafold_energy(utr3_sequence, rnafold_command)
                    txfeat_dict.update({'utr3.structure_mfe': region_energy['mfe'],
                                        'utr3.structure_mea': region_energy['mea'],
                                        'utr3.structure_centroid': region_energy['centroid']})
        except:
            logger.debug('Job %s Error: Transcript %s three_prime_UTR feature undetermined' % (job, tx_read.tx_id))

        # get selenocysteine feature
        try:
            if tx_read.tx_status['stop_codon_redefined_as_selenocysteine'] == 'defined':
                txfeat_dict.update({'stop_codon.seleno': 'yes'})
            else:
                txfeat_dict.update({'stop_codon.seleno': 'no'})
        except:
            logger.debug('Job %s Error: Transcript %s stop_codon_seleno feature undetermined' % (job, tx_read.tx_id))

        # append feature to return list
        tx_feature.append(txfeat_dict)

        # progress bar up increment
        # progress += 1
        stepper += 1

    # Completion time
    task_time = format(round((time.time() - initial_time) / 60, 2), '0.2f')
    # if job == 0:
    #    utils.progress_bar(1, 1, status='Complete!')
    if job == chunk_len - 1:
        utils.waiting_bar(1, message='Chunk %i / %i completed' % (job, chunk_len))
        logger.debug('Feature aggregation took an average of %s minutes per chunk' % (task_time))

    # Return tx_assembled
    return tx_feature
