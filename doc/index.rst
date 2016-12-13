=======================
Liknorm's documentation
=======================

Estimates

.. math::

  \int_{l}^r \text{ExpFam}(\theta, g(x)) \mathcal{N} (x | \mu, \sigma^2) \mathrm d x

via deterministic numerical integration.

:math:`g(\cdot)` is a link function and :math:`\theta` determines the
exponential-family distribution of interest:

- Bernoulli
- Binomial
- Poisson
- Exponential
- Gamma
- Geometric

---------------------
Functions description
---------------------

.. c:function:: LikNormMachine* liknorm_create_machine(int size)

  Create a Machine instance capable of doing numerical integration.

  :param int size: Number of integration points. ``500`` points should be
                   enough. ``300`` is usually fine too.
  :return: Machine instance to perform integration.
  :rtype: LikNormMachine*

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
