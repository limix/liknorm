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
- Bernoulli (Probit) :c:func:`liknorm_set_probit`
- Binomial :c:func:`liknorm_set_binomial`
- Poisson :c:func:`liknorm_set_poisson`
- Exponential :c:func:`liknorm_set_exponential`
- Gamma :c:func:`liknorm_set_gamma`
- Geometric :c:func:`liknorm_set_geometric`

-------
Install
-------

You can install it via conda

.. code-block:: bash

  conda install -c conda-forge liknorm

or by cloning this repository and building it

.. code-block:: bash

  git clone https://github.com/glimix/liknorm.git
  cd liknorm
  mkdir build
  cd build
  cmake ..
  make
  sudo make install

-------------
Usage example
-------------

Suppose you have the file

.. code-block:: c

  /* example.c */
  #include "liknorm.h"

  #include <stdio.h>

  int main()
  {
    double log_zeroth, mean, variance;
    double prior_var = 2.5;
    double prior_mean = -2.0;
    double nsuccesses = 2;
    double ntrials = 15;

    struct LikNormMachine *machine = liknorm_create_machine(500);

    liknorm_set_binomial(machine, nsuccesses, ntrials);
    liknorm_set_prior(machine, 1 / prior_var, prior_mean / prior_var);

    liknorm_integrate(machine, &log_zeroth, &mean, &variance);

    printf("%f\n", log_zeroth);
    printf("%f\n", mean);
    printf("%f\n", variance);

    liknorm_destroy_machine(machine);
  }

Compiling, linking, and running it via

.. code-block:: bash

  cc libliknorm.a example.c -o example
  ./example

should print::

  -2.049961
  -2.038184
  0.524308

---------------------
Functions description
---------------------

.. c:function:: struct LikNormMachine* liknorm_create_machine(int size)

  Create a Machine instance capable of doing numerical integration.

  :param int size: Number of integration points. ``500`` points should be
                   enough. ``300`` is usually fine too.
  :return: Machine instance to perform integration.
  :rtype: struct LikNormMachine*

.. c:function:: void liknorm_integrate(struct LikNormMachine *machine, double *log_zeroth, double *mean, double *variance)

  Perform numerical integration.

  :param struct LikNormMachine* machine: Machine to perform integration.
  :param double* log_zeroth: Zeroth moment.
  :param double* log_mean: First moment of the normalized distribution.
  :param double* log_variance: Variance of the normalized distribution.

.. c:function:: void liknorm_destroy_machine(struct LikNormMachine *machine)

  Destroy a Machine instance.

  :param struct LikNormMachine* machine: Machine to be destroyed. Always call it before
                                 exiting your program, otherwise it will
                                 leak memory.

.. c:function:: void liknorm_set_bernoulli(struct LikNormMachine *machine, double k)

  Set a Bernoulli likelihood.

  :param struct LikNormMachine* machine: Machine to perform integration.
  :param double k: ``0`` or ``1`` indicating a Bernoulli outcome.

.. c:function:: void liknorm_set_probit(struct LikNormMachine *machine, double k)

  Set a Bernoulli (Probit) likelihood.

  :param struct LikNormMachine* machine: Machine to perform integration.
  :param double k: ``0`` or ``1`` indicating a Bernoulli outcome.

.. c:function:: void liknorm_set_binomial(struct LikNormMachine *machine, double k, double n)

  Set a Binomial likelihood.

  :param struct LikNormMachine* machine: Machine to perform integration.
  :param double k: Number of successes.
  :param double n: Number of trials.

.. c:function:: void liknorm_set_poisson(struct LikNormMachine *machine, double k)

  Set a Poisson likelihood.

  :param struct LikNormMachine* machine: Machine to perform integration.
  :param double k: Number of successes.

.. c:function:: void liknorm_set_exponential(struct LikNormMachine *machine, double x)

  Set a Exponential likelihood.

  :param struct LikNormMachine* machine: Machine to perform integration.
  :param double x: Time span.

.. c:function:: void liknorm_set_gamma(struct LikNormMachine *machine, double x, double a)

  Set a Gamma likelihood.

  :param struct LikNormMachine* machine: Machine to perform integration.
  :param double x: Positive outcome.
  :param double a: Shape parameter.

.. c:function:: void liknorm_set_geometric(struct LikNormMachine *machine, double x)

  Set a Geometric likelihood.

  :param struct LikNormMachine* machine: Machine to perform integration.
  :param double x: Number of trials to success.

.. c:function:: void liknorm_set_prior(struct LikNormMachine *machine, double tau, double eta)

  Set the natural parameters of Normal prior.

  :param struct LikNormMachine* machine: Machine to perform integration.
  :param double tau: It equals to :math:`\sigma^{-2}`.
  :param double eta: It equals to :math:`\mu \sigma^{-2}`.
