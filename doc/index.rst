=======================
Liknorm's documentation
=======================

.. c:function:: LikNormMachine* liknorm_create_machine(int size)

  Create a Machine instance capable of doing numerical integration.

  :param int size: Number of integration points. ``500`` points should be
                   enough. ``300`` is usually fine too.
  :return: Machine instance to perform integration.
  :rtype: LikNormMachine

.. c:function:: void liknorm_destroy_machine(LikNormMachine *machine)

  Destroy a Machine instance.

  :param LikNormMachine machine: Machine to be destroyed. Always call it before
                                 exiting your program, otherwise it will
                                 leak memory.

.. c:function:: void liknorm_set_bernoulli(LikNormMachine *machine, double k)

  Set a Bernoulli likelihood.

  :param LikNormMachine machine: Machine to perform integration.
  :param double k: ``0`` or ``1`` indicating a Bernoulli outcome.

.. c:function:: void liknorm_set_prior(LikNormMachine *machine, double tau, double eta)

  Set the natural parameters of Normal prior.

  :param LikNormMachine machine: Machine to perform integration.
  :param double tau: It equals to :math:`\sigma^{-2}`.
  :param double eta: It equals to :math:`\mu / \sigma^2`.
