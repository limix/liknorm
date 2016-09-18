# liknorm

C library for computing the mean and variance of the product of an
exponential-family likelihood with a Normal distribution.

Given any exponential-family distribution $p(y|x)$, we want to compute

\[
\int_{-\infty}^{\infty} x^m p(y|x) \mathcal N(x | \mu, \sigma^2)
    \mathrm dx / Z,
\]

where $Z = \int_{-\infty}^{\infty} p(y|x) \mathcal N(x | \mu, \sigma^2) dx$
normalizes it and $m \in \{1, 2\}$, in less than a millisecond.
