"""
Checks system for required software needed to run the build_db pipeline.
"""
import sys
import logging
from shutil import which
from typing import List


def ready():
    """Check if software in list are in PATH and marked as executable."""

    # Setup logger
    logger = logging.getLogger(__name__)

    # Check if packages in list exists in path and executable
    packages = ['RNAfold', 'RNALfold']  # type: List[str]
    for wares in packages:
        if which(wares) is None:
            logger.info('Warning: %s could not be found! Please install before running build_db.' % wares)
            sys.exit()

    logger.debug('System requirements are satisfied, proceeding with building.')


if __name__ == '__main__':
    ready()
