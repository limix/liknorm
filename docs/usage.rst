*************
Usage example
*************

Suppose you have the file

.. code-block:: c

  /* example.c */
  #include "liknorm/liknorm.h"

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
