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
    keys = ['rnafold_command', 'rnalfold_command']
    for item in keys:
        if item not in list(cfg.keys()):
            return 'fail', item
        if item == 'viennarna_dir':
            if cfg[item][0] == '~':
                cfg[item] = cfg[item].replace('~', env_variables.base_path)
    return 'pass', 0


# Load and read configuration file
def load_check():
    # setup logging
    logger = logging.getLogger(__name__ + '.load_check')

    if os.path.isfile(env_variables.cfg_path + 'build_db.config'):
        cfg_file = open(env_variables.cfg_path + 'build_db.config', 'r')
        config_parse(cfg_file)
        status, key = config_check()
        if status == 'fail':
            logger.info('Error: missing %s in build_db.config file!' % key)
            sys.exit(1)
        else:
            logger.debug('Configurations from build_db.config loaded.')
    else:
        logger.info('Error: build_db.config file not found!')
        sys.exit(1)
