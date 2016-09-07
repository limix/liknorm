from numpy.testing import assert_almost_equal

from limix_math.special import logbinom

def test_lmath_logbinom():
    assert_almost_equal(logbinom(10, 3), 4.78749174278204669974)
