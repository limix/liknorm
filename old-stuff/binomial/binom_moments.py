from . import _binomial_ffi
from numba import cffi_support as _cffi_support
_cffi_support.register_module(_binomial_ffi)

from numpy import empty
from numpy import ndarray

def ptr(a):
    if isinstance(a, ndarray):
        return _binomial_ffi.ffi.cast("double *", a.ctypes.data)
    return a

class BinomialMoments(object):
    def __init__(self, nintervals):
        super(BinomialMoments, self).__init__()
        self._nintervals = nintervals
        self._height = empty(nintervals+1)
        self._weight = empty(nintervals)

    def compute(self, nsuc, ntrials, eta, tau, lmom0, mu, var):
        ma = _binomial_ffi.lib.moments_array

        args = [nsuc, ntrials, eta, tau, lmom0, mu, var, len(nsuc),
                self._nintervals, self._height, self._weight]

        args = [ptr(a) for a in args]

        ma(*args)
