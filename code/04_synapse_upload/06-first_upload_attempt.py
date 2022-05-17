#  This script is designed to perform an upload to synapse in a single command,
#  using the synapse python client. A user passes the path to the "manifest"
#  and a path to a YAML file containing synapse username, password, and
#  personal access token. To perform a dry run, pass the "-d" argument.
#
#  This script can be invoked like this:
#
#      python 05_first_upload_attempt.py \
#          -c <path to YAML file containing synapse credentials> \
#          -m <path to manifest> \
#          -d <dry run? (default: False)>
#
#  In practice, the synapse/2.6.0 module is loaded to set the version of python
#  and required dependencies.
#
#  ---
#  username: "some_user"
#  password: "some_password"
#  token: "some_token"
#

import sys
import getopt
import synapseclient
import synapseutils
import yaml

#  Recieve command-line options
try:
    opts, args = getopt.getopt(
        sys.argv[1:], "c:m:d", ["credentials=", "man_path=", "dry_run"]
    )
except getopt.GetoptError:
    print('python 05_first_upload_attempt.py \\')
    print('    -c <path to YAML file containing synapse credentials> \\')
    print('    -m <path to manifest>')
    sys.exit(2)

dry_run = False

for opt, arg in opts:
    if opt in ('-c', '--credentials'):
        #  Read in YAML file of credentials
        with open(arg, 'r') as f:
            cred_dict = yaml.safe_load(f)
        
        #  Grab values from YAML file
        for key, value in cred_dict.items():
            assert key in ('username', 'password', 'token'), \
                key + ' is not a supported key. Only "username", "password", and "token" are accepted.'
                
            if key == 'username':
                username = value
            elif key == 'password':
                synapse_password = value
            else: # must be 'token'
                token = value
    
    elif opt in ('-m', '--man_path='):
        man_path = arg
    elif opt in ('-d', '--dry_run'):
        dry_run = True
    else:
        sys.exit('Unhandled option.')


#  Log in to synapse
syn = synapseclient.Synapse()
syn.login(authToken=token)

if dry_run:
    print('Performing a dry run only!')

#  Upload
synapseutils.sync.syncToSynapse(syn, man_path, dryRun=dry_run)
