Estimates

.. math::

    \int_{-\inf}^{\inf} \text{ExpFam}(\theta, g(x)) \mathcal{N} (x | \mu, \sigma^2) \mathrm d x

via a very fast and deterministic numerical integration.

*******
Install
*******

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

.. toctree::
    :caption: Table of contents
    :name: mastertoc
    :maxdepth: 5

    introduction
    usage
    interface
    implementation

*****************
Comments and bugs
*****************

You can get the source and open issues `on Github.`_

.. _on Github.: https://github.com/limix/liknorm
