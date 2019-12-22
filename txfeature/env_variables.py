"""
contains environment variables
"""

import os

#base_path = "".join(str(os.path.dirname(os.path.realpath(__file__))).split('txfeature/')[0:2])
base_path = str(os.path.dirname(os.path.realpath(__file__)))
db_builder_path = base_path + "/db_builder/"
bin_path = base_path + "/bin/"
test_path = base_path + "/tests/"
