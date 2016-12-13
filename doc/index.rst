=======================
Liknorm's documentation
=======================

.. c:function:: LikNormMachine* liknorm_create_machine(int size)

  :param int size: Number of integration points. `500` points should be
                   enough.
  :return: Machine instance to perform integration.
  :rtype: LikNormMachine

.. c:function:: void liknorm_destroy_machine(LikNormMachine *machine)

  :param LikNormMachine machine: Machine to be destroyed. Always call it before
                                 exiting your program, otherwise it will
                                 leak memory.

.. c:function:: void liknorm_set_bernoulli(LikNormMachine *machine, double k)

  :param LikNormMachine machine: Machine to perform integration.
  :param double k: `0` or `1` indicating a Bernoulli outcome.
