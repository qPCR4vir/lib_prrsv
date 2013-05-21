#
# Nagios interface data and functions
#
__version__ = [int(x) for x in "$Revision: 1.2 $".split()[1].split('.')]
__author__ = "Kael Fischer <kael.fischer@gmail.com>"

import getopt as G
import sys

OK = 0
WARN = 1
CRITICAL = 2

state={0:'OK',
       1:'WARN',
       2:'CRITICAL'
       }

def print4nagios (status,message):
    """Produce and output the nagios oneliner
    """
    print "%s: %s" %(state[status],message)

def getopt(scriptArgs,pluginDoc):
    """process nagions plugin opts
    """
    try:
        optlist, args = G.getopt(scriptArgs,'hH?')
    except:
        raise
        print
        print usage(pluginDoc)
        sys.exit(1)
        
    for (o,a) in optlist:
        if o in ('-h','-H','-?'):
            print
            print usage(pluginDoc)
            sys.exit(0)
            
    if len(args) != 1:
        print
        print usage(pluginDoc)
        sys.exit(2)
    else:
        return args[0]

def usage(pluginDoc):
    """the default nagios plug-in usage message
    """

    return ("""USAGE %s [-Hh]|host_IP

OPTIONS
    -h    Print this message
    -H    Print this message

PLUGIN DOCUMENTATION
""" %(sys.argv[0])) + pluginDoc
            

    
