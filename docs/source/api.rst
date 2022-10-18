API
===
    
Abaqus Interface
----------------

.. doxygenfunction:: umat
   :project: umat

Abstract Classes
----------------

The following classes form the basis for user-implemented models.

Model 
^^^^^

The base class from which all constitutive models are derived.

.. doxygenclass:: Model
    :project: umat
    :members:
    :protected-members:
    :private-members:

Elastic
^^^^^^^

Elastic model class that contains basic elastic functionality. Inherits the Model class.

.. doxygenclass:: Elastic
    :project: umat
    :members:
    :protected-members:
    :private-members:

Elastoplastic
^^^^^^^^^^^^^

Elastoplastic model class that contains basic elastoplastic functionality based on refined explicit stress integration with 
automatic error control after Sloan et al. (2001). Inherits the Model and Elastic classes.

.. doxygenclass:: Elastoplastic
    :project: umat
    :members:
    :protected-members:
    :private-members:

Models
------

The following models have been implemented using this framework. Each model inherits one of the "Abstract Classes" based 
on the genera of constitutive behaviour required.

Modified Cam Clay (MCC)
^^^^^^^^^^^^^^^^^^^^^^^

The Modified Cam Clay (MCC) model is implemented here. Theory to go here...

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

Utilities
---------

The following typdefs and classes contain types used by the constitutive model development framework. 

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