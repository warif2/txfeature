"""

"""

import datetime
import time
import argparse
import logging
import sys
import os
import csv
import concurrent.futures
import pandas as pd

from txfeature.db_builder import gff_parser, tx_assembly, tx_features, txfeat_functions, build_config
from txfeature.db_builder import system_check
from txfeature import version


def main():
    # Start time of analysis
    time_stamp = str(datetime.datetime.now())
    initial_time = time.time()

    # Setup of argparse for script arguments
    parser = argparse.ArgumentParser(
        description="Obtain a table of features associated to annotated transcripts within a gff3 annotation file.",
        prog="txfeature build_db")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-gff", type=str, default=None, metavar="GFF",
                          help="specify path to the associated gff3 file", required=True)
    required.add_argument("-fa", type=str, default=None, metavar="FASTA", help="specify path to the fasta file",
                          required=True)
    required.add_argument("-out", type=str, default=None, metavar="OUTPUT", help="label for output directory",
                          required=True)
    optional.add_argument("-t", "--threads", nargs='?', const=1, type=int, default=1,  metavar="",
                          help='number of threads to utilize (default = 1)')
    optional.add_argument("-v", "--version", action='version', version='%(prog)s v' + version.__version__)
    optional.add_argument("-s", action='store_true', help="silence terminal output")
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # Creating sub-directories in output path
    sub_dir = ['log', 'txfeat_db']
    if not os.path.exists(args.out):
        os.makedirs(args.out)
    for folder in sub_dir:
        if not os.path.exists(args.out + '/' + folder):
            os.makedirs(args.out + '/' + folder)

    # Preparing logging console for __main__
    numeric_level = getattr(logging, 'INFO', None)
    logging.basicConfig(filename=args.out + '/log/txfeature.' + time_stamp.replace(" ", "_") + '.log',
                        level=logging.DEBUG,
                        format='%(asctime)s\t%(name)-12s\t%(message)s',
                        filemode='w')
    logger = logging.getLogger('db_builder.main')
    logger.debug('db_builder.py version: %s' % version.__version__)
    logger.debug('Input command: python db_builder.py ' + " ".join(sys.argv))

    # Defining Handler to write messages to sys.stdout
    if args.s:
        sys.stdout = open(os.devnull, 'w')
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(numeric_level)
    formatter = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%y-%m-%d %H:%M:%S')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)
    logger.info('Starting database construction with txfeature ver=%s' % version.__version__)

    # Load and read configuration file
    logger.debug('Performing system check to ensure necessary executables are installed.')
    system_check.ready()
    build_config.load_check()

    # Parse gff into searchable dataframe
    logger.info('Parsing gene annotation file...')
    gff_num = sum(1 for _ in open(args.gff))
    gff_df = gff_parser.gff_table(csv.reader(open(args.gff, 'r'), delimiter='\t'))
    logger.info('Parsing complete!')
    logger.info('Number of entries in gff: %i' % gff_num)
    logger.info('Number of annotated transcripts: %i' % len(gff_df['tx_attr'].keys()))

    # Setup parallel processing
    logger.info('Preparing transcript assembly for %i threads.' % args.threads)
    chunk_size = int(len(gff_df['tx_attr'])/args.threads)
    tx_chunks = txfeat_functions.chunks(gff_df['tx_attr'], chunk_size)
    tx_jobs = []
    for item in tx_chunks:
        tx_jobs.append(item)

    # Start tx assembly jobs
    tx_assembled = {}
    logger.info('Starting assembly of transcripts...')
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        njobs = range(len(tx_jobs))
        jobs = [executor.submit(tx_assembly.assemble, gff_df, tx_jobs[i], args.fa, args.out, i) for i in njobs]

        # collect results from jobs
        for job in concurrent.futures.as_completed(jobs):
            tx_assembled.update(job.result())
    logger.debug('Transcript assembly complete!')

    # Construct txfeat dataframe
    txfeat_df = []
    logger.info('Aggregating features for transcripts...')
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        jobs = [executor.submit(tx_features.add_features, tx_assembled, tx_jobs[i], i) for i in njobs]

        # collect results from jobs
        for job in concurrent.futures.as_completed(jobs):
            txfeat_df += job.result()
    logger.debug('Transcript feature aggregation complete!')

    # Save results to csv using pandas
    logger.info('Saving txfeat_db to output directory...')
    output = pd.DataFrame(txfeat_df).fillna(pd.np.nan)
    output.to_csv(args.out + '/txfeat_db_' + args.out + '.csv', index=False, na_rep='NULL')

    # Completion time
    task_time = format(round((time.time() - initial_time) / 60, 2), '0.2f')
    logger.debug('It took %s minutes to construct txfeature database' % task_time)


if __name__ == '__main__':
    main()
