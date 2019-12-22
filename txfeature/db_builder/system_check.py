"""
Checks system for required software needed to run the build_db pipeline.
"""
import sys
import os
import subprocess
import logging
from shutil import which
from typing import List
from txfeature.db_builder import build_config


def ready():
    """Check if software in list are in PATH and marked as executable."""

    # Setup logger
    logger = logging.getLogger(__name__)

    # Check if packages in list exists in path and executable
    packages = ['RNAfold', 'RNALfold']  # type: List[str]
    if build_config.cfg['viennarna_dir'] == '':
        for wares in packages:
            if which(wares) is None:
                logger.info('Warning: %s could not be found! Please install before running build_db.' % wares)
                sys.exit()
    else:
        if os.path.isdir(build_config.cfg['viennarna_dir']):
            for wares in packages:
                try:
                    cmd = build_config.cfg['viennarna_dir'] + wares + ' -h'
                    subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
                except subprocess.CalledProcessError:
                    logger.info('Warning: %s could not be found! Please install before running build_db.' % wares)
                    sys.exit()

    logger.debug('System requirements are satisfied, proceeding with building.')


if __name__ == '__main__':
    ready()
