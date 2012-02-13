import time
from datetime import timedelta
from functools import wraps

from credovi import app

class Timer(dict):
    """
    Timer used to benchmark sections of code in application debugging mode.
    """
    def start(self, name='default'):
        """
        """
        self.__setitem__(name, time.clock())
    
    def elapsed(self, name='default'):
        """
        """
        return time.clock() - self.__getitem__(name)
    
    def formatted(self, name='default'):
        """
        """
        return timedelta(seconds=self.elapsed(name))

def timer(message):
    '''
    Decorator that is used to log the execution time of a function.
    '''
    def wrap(function):
        def wrapped(*args, **kwargs):
            elapsed = time.clock()
            result = function(*args, **kwargs)
            app.log.debug(message.format(time.clock() - elapsed))
            
            return result

        return wraps(function)(wrapped)
    return wrap