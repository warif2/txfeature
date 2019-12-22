"""

"""

# Modules
import sys
import shutil
import db_builder.db_builder as db_build
import txfeature.env_variables

help_text_txfeat = """txfeature v1.0
A tool to obtain features associated to annotated transcripts within
gff files.

usage: txfeature <mode> [options] [-h]

modes:
build_db          pipeline to construct transcript feature database
build_db_config   output configuration file to working directory to modify settings

***Further [options] for mode can be obtained by <mode> --help."""

if __name__ == '__main__':

    # Run db_builder.py
    if sys.argv[1] == 'build_db':
        sys.argv = sys.argv[1:]
        db_build.main()

    elif sys.argv[1] == 'build_db_config':
        cfg = txfeature.env_variables.db_builder_path + 'build_db.cfg'
        shutil.copyfile(cfg, './build_db.cfg')

    elif sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print(help_text_txfeat)

    else:
        print('usage: txfeature <mode> [options] [-h] ')
        print("txfeature: error: unrecognized command '%s'" % sys.argv[1])
        print('Please refer to --help for usage and options.')
