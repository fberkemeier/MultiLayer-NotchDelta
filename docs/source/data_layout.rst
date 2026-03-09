Data Layout
===========

Expected structure
------------------

Provide data under ``data/`` using this layout:

- ``data/adjacency_matrices/adjacency_matrices_<wing_region>.xlsx``
- ``data/cell_geometry/cell_geometry_<wing_region>.csv``
- ``data/wing_region_metadata.csv``

Naming rule
-----------

For each metadata ``wing_region`` value, matching adjacency and geometry files must exist.
For example, ``wing_region = wd_1`` maps to:

- ``data/adjacency_matrices/adjacency_matrices_wd_1.xlsx``
- ``data/cell_geometry/cell_geometry_wd_1.csv``

Adjacency matrices (3D topology input)
--------------------------------------

Each ``adjacency_matrices_<wing_region>.xlsx`` file contains one sheet per depth layer.
Each sheet stores a square adjacency matrix:

- rows and columns are cell indices (same ordering as model indices)
- value ``1`` means the cell pair shares a contact in that layer
- value ``0`` means no contact
- diagonal entries should be ``0`` (no self-contact)

The number of sheets should match ``n_layers`` in metadata for that region.

Cell geometry file
------------------

Each ``cell_geometry_<wing_region>.csv`` file stores cell-level geometry used by the model
(for centroids, apical area, and derived diameters).

Observed columns in current dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Example ``cell_geometry`` columns and meanings
   :header-rows: 1
   :widths: 24 56 20

   * - Column
     - Meaning
     - Example
   * - ``frame``
     - Layer/time identifier from source segmentation (typically ``0`` for static inputs)
     - ``0``
   * - ``label``
     - Cell index used in modelling arrays
     - ``4``
   * - ``area``
     - Cell area (apical cross-section units from source data)
     - ``1148.0``
   * - ``centroid_x``
     - X coordinate of cell centroid
     - ``172.639``
   * - ``centroid_y``
     - Y coordinate of cell centroid
     - ``12.360``
   * - ``original_label``
     - Original segmentation label before remapping
     - ``5.0``
   * - ``perimeter``
     - Cell perimeter from segmentation
     - ``85.0``
   * - ``show``
     - Visibility/filter flag from source preprocessing
     - ``1.0``

Wing region metadata
--------------------

``data/wing_region_metadata.csv`` defines region-level model parameters.

Required columns
^^^^^^^^^^^^^^^^

.. list-table:: Metadata columns and meaning
   :header-rows: 1
   :widths: 24 56 20

   * - Column
     - Meaning
     - Example
   * - ``wing_region``
     - Region key used in filenames and model calls
     - ``wd_1``
   * - ``wing_label``
     - Display label used in figures/legends
     - ``1``
   * - ``gap``
     - Inter-layer spacing (``?L``)
     - ``0.5``
   * - ``n_layers``
     - Number of depth layers used in simulation
     - ``50``
   * - ``default_height``
     - Default SOP-band height for spacing analyses
     - ``85``
   * - ``is_wing_disc``
     - Boolean selector for canonical wing-disc entries
     - ``TRUE``
   * - ``signalling_labels``
     - Serialized list of signalling-capable cell indices
     - ``[23, 25, 26, ...]``

Notes on ``signalling_labels``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Store ``signalling_labels`` as a JSON-like list string in each metadata row,
for example ``"[1, 5, 9, 12]"``.

