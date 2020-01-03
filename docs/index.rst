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


or by building it yourself via the following command:

.. code-block:: bash

    # DO_CMD=sudo
    curl -fsSL https://git.io/JerYI | GITHUB_USER=horta GITHUB_PROJECT=logaddexp bash

The above command should work on Windows, Linux, and MacOS.

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
