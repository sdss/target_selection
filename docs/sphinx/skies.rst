.. _skies:

Sky catalogues
==============

This section describes how parent sky catalogues are created and how they are used to generate sky cartons.


Creating the parent sky catalogue
---------------------------------

Parent sky catalogues are generated as a list of sky position that are considered blank enough to be used for sky subtraction, i.e., that they are far enough from any neighbouring bright source. For a given source catalogue (e.g., Gaia) these steps are followed to generate a list of valid sky position (for that particular catalogue):

- The whole sky is divided in "tiles" of a given HEALPix order (usually :math:`k=5`). This is mostly done for the purposes of parallelising the code, so that each tile can be processed independently. For
- Each tile is further subdivided in
