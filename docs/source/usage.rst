Usage
=====

Python
------

The Python bindings allow rapid prototyping of single-element behaviour. After building and installing
the Python modules:

.. code-block:: python

    from geomat.models import MCC
    from geomat.utilities import Settings
    import numpy as np

    # Define parameters: M, nu, N, lambda*, kappa*
    parameters = np.array([0.92, 0.2, 1.195, 0.08, 0.02])
    state = np.array([10.0])  # Initial p_c

    # Create model.
    model = MCC(parameters, state)

    # Set initial stress state (compression-positive).
    model.set_sigma_prime_tilde(np.array([7.5, 7.5, 7.5, 0.0, 0.0, 0.0]))

    # Apply strain increment.
    model.set_Delta_epsilon_tilde(np.array([0.0005, -0.00025, -0.00025, 0.0, 0.0, 0.0]))
    model.solve()

    # Retrieve results.
    sigma = model.get_sigma_prime_tilde()
    p = model.get_p_prime()
    q = model.get_q()

Abaqus
------

The geomat framework produces a precompiled shared library ``libstandardU.so`` (Linux) or
``standardU.dll`` (Windows) that can be used directly with Abaqus without runtime compilation.

Building
^^^^^^^^

.. code-block:: console

    mkdir build && cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make standardU

This produces ``libstandardU.so`` in the build directory.

Configuration
^^^^^^^^^^^^^

To use the precompiled library, set the ``usub_lib_dir`` variable in your Abaqus environment file
(``abaqus_v6.env``) to point to the directory containing the shared library:

.. code-block:: python

    usub_lib_dir = "/path/to/geomat/build"

Material Definition
^^^^^^^^^^^^^^^^^^^

In your Abaqus input file or CAE, define a user material with the appropriate model name
(e.g. ``MCC``, ``C2MC``, ``LinearElastic``, ``EMC``, ``SMCC``) and provide the parameters
in the ``*USER MATERIAL`` card. For example, for the Modified Cam Clay model:

.. code-block:: text

    *MATERIAL, NAME=MCC
    *USER MATERIAL, CONSTANTS=5, UNSYMM
    0.92, 0.2, 1.195, 0.08, 0.02

Note that Abaqus uses a tension-positive sign convention; the interface handles the conversion
to the compression-positive convention used internally.

Plaxis
------

The framework produces a shared library ``libusrmod.so`` (Linux) or ``usrmod.dll`` (Windows)
compatible with the Plaxis User-Defined Soil Model (UDSM) interface.

Building
^^^^^^^^

.. code-block:: console

    mkdir build && cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make usrmod

On Windows, compile to produce ``usrmod.dll``.

Installation
^^^^^^^^^^^^

Place the compiled DLL in the ``udsm`` sub-folder of the Plaxis program directory. Plaxis will
automatically detect the available models and their parameters via the metadata functions exported
by the DLL.

Material Selection
^^^^^^^^^^^^^^^^^^

In Plaxis Input, select *User-defined* from the Material model dropdown in the General tabsheet.
The Parameters tabsheet will show the available DLLs and models. Select the appropriate model and
enter the parameter values. Parameter names and units are read automatically from the DLL.

Note that Plaxis uses a compression-negative sign convention; the interface handles the conversion
to the compression-positive convention used internally.
