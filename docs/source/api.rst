API
===
    
Abaqus Interface
----------------

.. doxygenfunction:: umat
   :project: geomat

Abstract Classes
----------------

The following classes form the basis for user-implemented models.

Model 
^^^^^

The base class from which all constitutive models are derived.

.. doxygenclass:: Model
    :project: geomat
    :members:
    :protected-members:
    :private-members:

Elastic
^^^^^^^

Elastic model class that contains basic elastic functionality. Inherits the Model class.

.. doxygenclass:: Elastic
    :project: geomat
    :members:
    :protected-members:
    :private-members:

Elastoplastic
^^^^^^^^^^^^^

Elastoplastic model class that contains basic elastoplastic functionality based on refined explicit stress integration with 
automatic error control after Sloan et al. (2001). Inherits the Model and Elastic classes.

.. doxygenclass:: Elastoplastic
    :project: geomat
    :members:
    :protected-members:
    :private-members:

Models
------

The following models have been implemented using this framework. Each model inherits one of the "Abstract Classes" based 
on the genera of constitutive behaviour required.

Linear Elastic
^^^^^^^^^^^^^^

The Linear Elastic model is implemented here. Theory to go here...

.. doxygenclass:: LinearElastic
    :project: geomat
    :members:
    :protected-members:
    :private-members:

Modified Cam Clay (MCC)
^^^^^^^^^^^^^^^^^^^^^^^

The Modified Cam Clay (MCC) model is implemented here. Theory to go here...

.. doxygenclass:: MCC
    :project: geomat
    :members:
    :protected-members:
    :private-members:

Soft Modified Cam Clay (SMCC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Soft Modified Cam Clay (SMCC) model is implemented here. Theory to go here...

.. doxygenclass:: SMCC
    :project: geomat
    :members:
    :protected-members:
    :private-members:

Utilities
---------

The following modules, typedefs and classes contain types used by the constitutive model development framework. 

Checks
^^^^^^

.. doxygennamespace:: Checks
    :project: geomat
    :members:

Logging
^^^^^^^

.. doxygennamespace:: Logging
    :project: geomat
    :members:

Intersection
^^^^^^^^^^^^

.. doxygennamespace:: Intersection
    :project: geomat
    :members:

Types
^^^^^

.. doxygentypedef:: Voigt
    :project: geomat

.. doxygentypedef:: Cauchy
    :project: geomat

.. doxygentypedef:: Constitutive
    :project: geomat

.. doxygentypedef:: Jacobian
    :project: geomat

.. doxygentypedef:: Parameters
    :project: geomat

.. doxygentypedef:: State
    :project: geomat

.. doxygenfunction:: to_cauchy
    :project: geomat

.. doxygenfunction:: to_voigt
    :project: geomat