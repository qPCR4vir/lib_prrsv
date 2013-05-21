#!/usr/local/bin/python
"""utilities for turning module __doc__'s into
optparse.OptionParsers.

adapted from: M. Simionato, Henry Crutcher
              (http://code.activestate.com/recipes/278844/)

"""
__verison__ = tuple([int(x) for x in
                     '$Revision: 1.2 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

import optparse
import re
import sys

class ParsingError(Exception):
    pass

def moduleOptionParser(module):
    """Factory function to return an optparse.OptionParser istance based on the
    __doc__ attribute of the module.  Usage message also includes __version__
    and __author__ if defined.
    """
    USAGE = re.compile(r'(?s)\s*[Uu][Ss][Aa][Gg][eE]: (.*?)(\n[ \t]*\n|$)')

    docStr = module.__doc__
    try:
        author = "Author: %s" % module.__author__
    except:
        author=''
    try:
        version = 'Version: %s' % module.__version__
    except:
        version=''

    uMatch = USAGE.search(docStr)
    if uMmatch == None:
        raise ParsingError("Cannot find the option string")

    optLines = match.group(1).splitlines()

    rv = optparse.OptionParser('\n'.join((docStr,version,author)),add_help_option=False)
    try:
        for line in optlines[1:]:
            opt, help=line.split(':')[:2]
            short,long=opt.split(',')[:2]
            if '=' in opt:
                # unless the value on the other side of = is the same
                # (modulo case) it is used as the default

                action='store'
                long, default=long.split('=')[:2]
                if default.lower()==long:
                  default=None
            else:
                default=
                action='store_true'
            rv.add_option(short.strip(),long.strip(),
                          action = action, help = help.strip(), default=default)
    except (IndexError,ValueError):
        raise ParsingError("Cannot parse the option string correctly:\n%s" % line)
    return p.parse_args(arglist)



if __name__ == "__main__":
    print sys.modules['__main__'].__author__
    
