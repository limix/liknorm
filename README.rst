
liknorm
=======

|Build-Status| |Win-Build-Status| |Codacy-Grade| |Doc-Status|

C library for computing moments of the product of an
exponential-family likelihood with a Normal distribution.

Install
-------

You can install it via conda_

.. code:: bash

    conda install -c conda-forge liknorm

A second installation option would be to download the latest source and to
build it by yourself.
On Linux or macOS systems it can be as simple as

.. code:: bash

    bash <(curl -fsSL https://raw.githubusercontent.com/limix/liknorm/master/install)

Refer to documentation_ for more information.

Authors
-------

* `Danilo Horta`_

License
-------

This project is licensed under the MIT License - see the `License file`_ file
for details.


.. |Build-Status| image:: https://travis-ci.org/limix/liknorm.svg?branch=master
    :target: https://travis-ci.org/limix/liknorm

.. |Win-Build-Status| image:: https://ci.appveyor.com/api/projects/status/kb4b4rcsm4t60bg5/branch/master?svg=true
    :target: https://ci.appveyor.com/project/Horta/liknorm/branch/master

.. |Codacy-Grade| image:: https://api.codacy.com/project/badge/Grade/689b555393364226863c3a237f801650
    :target: https://www.codacy.com/app/danilo.horta/liknorm?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=limix/liknorm&amp;utm_campaign=Badge_Grade

.. |Doc-Status| image:: https://readthedocs.org/projects/liknorm/badge/?style=flat-square&version=stable
    :target: https://liknorm.readthedocs.io/

.. _conda: http://conda.pydata.org/docs/index.html

.. _License file: https://raw.githubusercontent.com/limix/liknorm/master/LICENSE.txt

.. _Danilo Horta: https://github.com/horta

.. _documentation: http://liknorm.readthedocs.io/
