"""
Load and check build_db.config file.
"""

import os
from txfeature import env_variables
import logging
import sys

cfg = {}


def config_parse(file):
    global cfg
    for line in file:
        if line[0] == "#":
            continue
        else:
            configs = line.strip('\n').split('=')
            if len(configs) == 2:
                cfg[configs[0]] = configs[1]


def config_check():
    global cfg
    keys = ['rnafold_command', 'rnalfold_command', 'viennarna_dir']
    for item in keys:
        if item not in list(cfg.keys()):
            return 'fail', item
    return 'pass', 0


# Load and read configuration file
def load_check(file):
    # setup logging
    logger = logging.getLogger(__name__ + '.load_check')

    if file is None:
        if os.path.isfile(env_variables.db_builder_path + 'build_db.cfg'):
            cfg_file = open(env_variables.db_builder_path + 'build_db.cfg', 'r')
            config_parse(cfg_file)
            status, key = config_check()
            if status == 'fail':
                logger.info('Error: missing %s in build_db.cfg file!' % key)
                sys.exit(1)
            else:
                logger.debug('Configurations from build_db.cfg loaded.')
        else:
            logger.info('Error: build_db.cfg file not found!')
            sys.exit(1)
    else:
        logger.debug('Loading user-defined configuration file.')
        if os.path.isfile(file):
            cfg_file = open(file, 'r')
            config_parse(cfg_file)
            status, key = config_check()
            if status == 'fail':
                logger.info('Error: missing %s in build_db.cfg file!' % key)
                sys.exit(1)
            else:
                logger.debug('Configurations from build_db.cfg loaded.')
        else:
            logger.info('Error: build_db.cfg file not found!')
            sys.exit(1)
