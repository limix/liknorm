=======================
Liknorm's documentation
=======================

Estimates

.. math::

  \int_{l}^r \text{ExpFam}(\theta, g(x)) \mathcal{N} (x | \mu, \sigma^2) \mathrm d x

via deterministic numerical integration.

The values :math:`l` and :math:`r` are necessarily the lower and upper bound
of the likelihood domain. (The user is not allowed to set it.)
:math:`g(\cdot)` is a link function and :math:`\theta` determines the
exponential-family distribution of interest:

- Bernoulli :c:func:`liknorm_set_bernoulli`
- Binomial :c:func:`liknorm_set_binomial`
- Poisson :c:func:`liknorm_set_poisson`
- Exponential :c:func:`liknorm_set_exponential`
- Gamma :c:func:`liknorm_set_gamma`
- Geometric :c:func:`liknorm_set_geometric`

----------------
How to build it?
----------------

In the project folder, type

.. code-block:: bash

  mkdir build
  cd build
  cmake ..
  make

It should create the library

.. code-block:: bash

  liknorm/libliknorm.[extension]

And if you want to test it:

.. code-block:: bash

  make test

-----
Usage
-----

.. code-block:: c

  #include "liknorm/liknorm.h"

  int main()
  {
    double log_zeroth, mean, variance;
    double prior_var = 2.5;
    double prior_mean = -2.0;
    double nsuccesses = 2;
    double ntrials = 15;

    LikNormMachine *machine = liknorm_create_machine(500);

    liknorm_set_binomial(machine, nsuccesses, ntrials);
    liknorm_set_prior(machine, 1 / prior_var, prior_mean / prior_var);

    liknorm_integrate(machine, &log_zeroth, &mean, &variance);

    liknorm_destroy_machine(machine);
  }

---------------------
Functions description
---------------------

.. c:function:: LikNormMachine* liknorm_create_machine(int size)

  Create a Machine instance capable of doing numerical integration.

  :param int size: Number of integration points. ``500`` points should be
                   enough. ``300`` is usually fine too.
  :return: Machine instance to perform integration.
  :rtype: LikNormMachine*

.. c:function:: void liknorm_integrate(LikNormMachine *machine, double *log_zeroth, double *mean, double *variance)

  Perform numerical integration.

  :param LikNormMachine* machine: Machine to perform integration.
  :param double* log_zeroth: Zeroth moment.
  :param double* log_mean: First moment of the normalized distribution.
  :param double* log_variance: Variance of the normalized distribution.

.. c:function:: void liknorm_destroy_machine(LikNormMachine *machine)

  Destroy a Machine instance.

  :param LikNormMachine* machine: Machine to be destroyed. Always call it before
                                 exiting your program, otherwise it will
                                 leak memory.

.. c:function:: void liknorm_set_bernoulli(LikNormMachine *machine, double k)

  Set a Bernoulli likelihood.

  :param LikNormMachine* machine: Machine to perform integration.
  :param double k: ``0`` or ``1`` indicating a Bernoulli outcome.

.. c:function:: void liknorm_set_binomial(LikNormMachine *machine, double k, double n)

  Set a Binomial likelihood.

  :param LikNormMachine* machine: Machine to perform integration.
  :param double k: Number of successes.
  :param double n: Number of trials.

.. c:function:: void liknorm_set_poisson(LikNormMachine *machine, double k)

  Set a Poisson likelihood.

  :param LikNormMachine* machine: Machine to perform integration.
  :param double k: Number of successes.

.. c:function:: void liknorm_set_exponential(LikNormMachine *machine, double x)

  Set a Exponential likelihood.

  :param LikNormMachine* machine: Machine to perform integration.
  :param double x: Time span.

.. c:function:: void liknorm_set_gamma(LikNormMachine *machine, double x, double a)

  Set a Gamma likelihood.

  :param LikNormMachine* machine: Machine to perform integration.
  :param double x: Positive outcome.
  :param double a: Shape parameter.

.. c:function:: void liknorm_set_geometric(LikNormMachine *machine, double x)

  Set a Geometric likelihood.

  :param LikNormMachine* machine: Machine to perform integration.
  :param double x: Number of trials to success.

.. c:function:: void liknorm_set_prior(LikNormMachine *machine, double tau, double eta)

  Set the natural parameters of Normal prior.

  :param LikNormMachine* machine: Machine to perform integration.
  :param double tau: It equals to :math:`\sigma^{-2}`.
  :param double eta: It equals to :math:`\mu \sigma^{-2}`.
