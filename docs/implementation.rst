**************
Implementation
**************

ExpFam
======

.. doxygenstruct:: ExpFam
   :members:

Likelihood
==========

We assume the canonical link function for every likelihood.

Bernoulli
---------

``y`` assumes ``1`` or ``0`` for failure.
We make use of the Binomial implementation. So, please, refer to the next section for
details.

.. doxygenfunction:: bernoulli_log_partition
.. doxygenfunction:: bernoulli_log_partition_fderivative
.. doxygenfunction:: bernoulli_log_partition_derivatives

Binomial
--------

The random variable is given by ``y = k/n``. The support is therefore
``y ϵ {0/n, 1/n, ..., r/n}``. The exponential family functions are::

    𝜙      = n
    a(𝜙)   = 1/𝜙
    b(θ)   = log(1 + exp(θ))
    c(y,𝜙) = log(binom(n, y𝜙))

Let us define::

    𝜇 = E[y] = p.

The canonical link function and its inverse are given by::

    canonical(𝜇)     = log(𝜇/(1+𝜇)) = η
    canonical_inv(η) = 1/(1 + exp(-η))

.. doxygenfunction:: binomial_log_partition
.. doxygenfunction:: binomial_log_partition_fderivative
.. doxygenfunction:: binomial_log_partition_derivatives

Negative Binomial
-----------------

The random variable is given by ``y = k/r``. The support is therefore
``y ϵ {0/r, 1/r, ..., r/r}``. The exponential family functions are::

    𝜙 = r
    a(𝜙) = 1/𝜙
    b(θ) = -log(1 - exp(θ))
    c(y,𝜙) = log(binom(y𝜙 + 𝜙 - 1, y𝜙))

Let us define::

    𝜇 = E[y] = p / (1 - p)

The canonical link function and its inverse are given by::

    canonical(𝜇)     = log(𝜇 / (1 + 𝜇)) = η
    canonical_inv(η) = exp(η) / (1 - exp(η))

.. doxygenfunction:: nbinomial_log_partition
.. doxygenfunction:: nbinomial_log_partition_fderivative
.. doxygenfunction:: nbinomial_log_partition_derivatives


Poisson
-------

The support is ``y ϵ {0, 1, ...}``. The exponential family functions are::

    𝜙      = 1
    a(𝜙)   = 𝜙
    b(𝜃)   = exp(𝜃)
    b'(𝜃)  = exp(𝜃)
    b'(𝜃)  = exp(𝜃)
    c(y,𝜙) = -log(y!)

Let us define::

    𝜇 = E[y] = λ,

for which ``λ`` is the Poisson distribution parameter. The canonical link function and
its inverse are given by::

    canonical(𝜇)     = log(𝜇 / (1 + 𝜇)) = η
    canonical_inv(η) = exp(η) / (1 - exp(η))

.. doxygenfunction:: poisson_log_partition
.. doxygenfunction:: poisson_log_partition_fderivative
.. doxygenfunction:: poisson_log_partition_derivatives
