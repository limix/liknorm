************
Introduction
************

The single-parameter exponential family is a class of distributions that can be
expressed as::

    f(y; θ, 𝜙) = exp{(yθ - b(θ))/a(𝜙) + c(y,𝜙)}.

for which ``𝜙`` is assumed to be known.
The definition of the functions ``a(.)``, ``b(.)``, and ``c(.)`` determines a
probabilistic distribution having the canonical parameter ``θ``. The expectation of
``y``, denoted here by ``𝜇``, determines the value of ``θ`` via the following relation::

    b'(θ) = 𝜇

Still, the value ``𝜇`` is often set indirectly via the natural parameter ``η``, which
relates to each other through a link function ``g(.)``::

    η = g(𝜇)

If ``g(.)`` is the so-called canonical function, we have the desirable equality::

    θ = η
