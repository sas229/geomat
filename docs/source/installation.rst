Installation
============

To install and build the umat development library simply:

.. code-block:: console

    git clone https://github.com/sas229/umat.git
    cd umat
    git submodule init src/external
    git submodule update --recursive
    mkdir build && cd build
    cmake ..
    make