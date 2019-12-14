"""
contains environment variables
"""

import os

base_path = "".join(str(os.path.dirname(os.path.realpath(__file__))).split('txfeature/')[0:2])
cfg_path = base_path + "/config/"
bin_path = base_path + "/bin/"
test_path = base_path + "/tests/"
lib_path = base_path + "/lib/"
