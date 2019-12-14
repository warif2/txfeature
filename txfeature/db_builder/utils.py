"""
Contains a set of functions that are used in various modules within the pipeline.
"""

import sys
import subprocess

def progress_bar(count, total, status=''):
    bar_len = 50
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    if status == 'Complete!':
        sys.stdout.write('[%s] %s%s ...%s' % (bar, percents, '%', status) + ' ' * 20 + '\n')
    else:
        sys.stdout.write('[%s] %s%s ...%s' % (bar, percents, '%', status) + ' ' * 20 + '\r')
        sys.stdout.flush()

def waiting_bar(stepper, message=''):
    symbol = "-\|/"
    step = (stepper % 4)
    sys.stdout.write('%s [%s] \r' % (message, symbol[step]))
    sys.stdout.flush()

def stdout_from_command(command):
    p = subprocess.Popen(command,
                         stdout = subprocess.PIPE,
                         shell = True)
    return iter(p.stdout.readline, b'')