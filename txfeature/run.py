"""

"""

# Modules
import sys
import db_builder.db_builder as db_build

if __name__ == '__main__':

    # Run db_builder.py
    if sys.argv[1] == 'build_db':
        sys.argv = sys.argv[1:]
        db_build.main()

    elif sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print('coming soon...')

    else:
        print('usage: txfeature <mode> [options] [-h] ')
        print("txfeature: error: unrecognized command '%s'" % sys.argv[1])
        print('Please refer to --help for usage and options.')
