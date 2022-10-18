API
===
    
Abaqus Interface
----------------

.. doxygenfunction:: umat
   :project: umat

Base Class
----------

.. doxygenclass:: Model
    :project: umat
    :members:
    :protected-members:
    :private-members:

Utilities
---------

The following typdefs and classes contain utility functions that can be used to develop constitutive models. 

Types
^^^^^

.. doxygentypedef:: Vector6d
   :project: umat

.. doxygentypedef:: Constitutive
   :project: umat

.. doxygenclass:: Cauchy
    :project: umat
    :members:
    :protected-members:
    :private-members:

.. doxygenclass:: Voigt
    :project: umat
    :members:
    :protected-members:
    :private-members:

Elastic
^^^^^^^

.. doxygenclass:: Elastic
    :project: umat
    :members:
    :protected-members:
    :private-members:

Model Classes
-------------

The following models have been implemented using this framework.

Modified Cam Clay (MCC)
^^^^^^^^^^^^^^^^^^^^^^^

The Modified Cam Clay (MCC) model is implemented here.

.. doxygenclass:: MCC
    :project: umat
    :members:
    :protected-members:
    :private-members:

.. Soft Modified Cam Clay (MCC)
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. The Soft Modified Cam Clay (MCC) model is implemented here.

.. .. doxygenclass:: SMCC
..     :project: umat
..     :members:
..     :protected-members:
..     :private-members:

