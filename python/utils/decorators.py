# initally from from http://wiki.python.org/moin/PythonDecoratorLibrary
__version__ = tuple([int(x) for x in
                      '$Revision: 1.8 $'.split()[1].split('.')])
__author__ = "Kael Fischer and Julio Carlos Menendez"

import warnings
import functools
import sys

def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used."""

    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.warn_explicit(
            "Call to deprecated function %(funcname)s." % {
                'funcname': func.__name__,
            },
            category=DeprecationWarning,
            filename=func.func_code.co_filename,
            lineno=func.func_code.co_firstlineno + 1
        )
        return func(*args, **kwargs)
    return new_func


## Usage examples ##
@deprecated
def my_func():
    pass

#@other_decorators_must_be_upper
@deprecated
def my_func():
    pass

class countcalls(object):
   "Decorator that keeps track of the number of times a function is called."

   __instances = {}

   def __init__(self, f):
      self.__f = f
      self.__numcalls = 0
      countcalls.__instances[f] = self

   def __call__(self, *args, **kwargs):
      self.__numcalls += 1
      return self.__f(*args, **kwargs)

   @staticmethod
   def count(f):
      "Return the number of times the function f was called."
      return countcalls.__instances[f].__numcalls

   @staticmethod
   def counts():
      "Return a dict of {function: # of calls} for all registered functions."
      return dict([(f, countcalls.count(f)) for f in countcalls.__instances])


import time
def benchmark(func):
    """
    A decorator that print (on stderr) the wallclock time
    function takes to execute.
    """
    import time
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        t = time.time()
        res = func(*args, **kwargs)
        print >>sys.stderr, "%s: %s sec"%(func.__name__, time.time()-t)
        return res
    return wrapper



# Retry decorator with exponential backoff
def retry(tries, delay=3, backoff=2):
  """Retries a function or method until it returns True.
  
  delay sets the initial delay, and backoff sets how much the delay should
  lengthen after each failure. backoff must be greater than 1, or else it
  isn't really a backoff. tries must be at least 0, and delay greater than
  0."""

  if backoff <= 1:
    raise ValueError("backoff must be greater than 1")

  tries = math.floor(tries)
  if tries < 0:
    raise ValueError("tries must be 0 or greater")

  if delay <= 0:
    raise ValueError("delay must be greater than 0")

  def deco_retry(f):
    def f_retry(*args, **kwargs):
      mtries, mdelay = tries, delay # make mutable

      rv = f(*args, **kwargs) # first attempt
      while mtries > 0:
        if rv == True: # Done on success
          return True

        mtries -= 1      # consume an attempt
        time.sleep(mdelay) # wait...
        mdelay *= backoff  # make future wait longer

        rv = f(*args, **kwargs) # Try again

      return False # Ran out of tries :-(

    return f_retry # true decorator -> decorated function
  return deco_retry  # @retry(arg[, ...]) -> true decorator


class memoized(object):
   """Decorator that caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned, and
   not re-evaluated.
   """
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      try:
         return self.cache[args]
      except KeyError:
         self.cache[args] = value = self.func(*args)
         return value
      except TypeError:
         # uncachable -- for instance, passing a list as an argument.
         # Better to not cache than to blow up entirely.
         return self.func(*args)
   def __repr__(self):
      """Return the function's docstring."""
      return self.func.__doc__


class LogPrinter(object):
    """LogPrinter class which serves to emulates a file object and logs
       whatever it gets sent to a Logger object at the INFO level."""
    def __init__(self):
        """Grabs the specific logger to use for logprinting."""
        self.ilogger = logging.getLogger('logprinter')
        il = self.ilogger
        logging.basicConfig()
        il.setLevel(logging.INFO)
    
    def write(self, text):
        """Logs written output to a specific logger"""
        self.ilogger.info(text)

def logprintinfo(func):
    """Wraps a method so that any calls made to print get logged instead"""
    @functools.wraps(func)
    def pwrapper(*arg):
        stdobak = sys.stdout
        lpinstance = LogPrinter()
        sys.stdout = lpinstance
        try:
            return func(*arg)
        finally:
            sys.stdout = stdobak
    return pwrapper
