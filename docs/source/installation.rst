Installation
============

To install and build this geotechnical constitutive model development library simply:

.. code-block:: console

    git clone https://github.com/sas229/geomat.git
    cd geomat
    git submodule init src/external
    git submodule update --init --recursive
    mkdir build && cd build
    cmake ..
    make