"""
Component of the db_builder pipeline, uses gff structure and iterates through and uses tx_build to build
transcripts.
"""

import logging
import time

from . import tx_build, txfeat_functions, utils


def assemble(gff_df, tx_dict, fasta, tmp_dir, job):
    """
    :param gff_df: parsed gff file using gff_parser
    :param tx_dict: list of transcript id for assembly
    :param fasta: fasta file for associated gff
    :param tmp_dir: temporary directory for transcript building
    :param job: integer value of the job
    """
    # setup logger and time
    logger = logging.getLogger(__name__ + '.assemble')
    logger.debug('Build job %i assembling annotated transcripts...' % job)
    initial_time = time.time()

    # setup return structure and tx list
    tx_assembled = {}
    tx_list = sorted(tx_dict.keys())

    # progress bar variables
    progress = 0
    total_tx = len(tx_list)

    # iterate through tx_list and construct table
    for tx in tx_list:
        # display progress of txfeat construction
        if job == 0:
            utils.progress_bar(progress, total_tx, status='txid: %s' % tx)

        # assemble tx
        tx_annot = txfeat_functions.tx2gff_lookup(gff_df, tx)
        tx_assembled[tx] = tx_build.build(tx, tx_annot, tx_dict[tx], fasta, tmp_dir + '/_tmp_txfeature')
        # progress bar up increment
        progress += 1

    # Completion time
    task_time = format(round((time.time() - initial_time) / 60, 1), '0.1g')
    if job == 0:
        utils.progress_bar(1, 1, status='Complete!')
    logger.debug('Build job %i took %s seconds to assemble %i annotated transcripts' % (job, task_time, total_tx))

    # Return tx_assembled
    return tx_assembled
