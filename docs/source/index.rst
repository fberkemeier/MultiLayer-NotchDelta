A Multilayer Model for Notch-Delta Signalling
======================================================================

The Multilayer Signalling Model (MSM) is a computational toolkit for simulating and analysing Notch-Delta signalling in three-dimensional epithelial tissues. It accounts for depth-resolved cell-cell contacts across apical and lateral surfaces, enabling systematic exploration of how tissue geometry and signalling range influence lateral inhibition and the emergence of spatial patterns of sensory organ precursor (SOP) cells.

.. toctree::
   :maxdepth: 1
   :caption: Contents

   overview
   quickstart
   data_layout
   function_reference
   release_notes

.. figure:: _static/img/msm_model.png
   :alt: Multilayer Signalling Model schematic placeholder
   :width: 80%
   :align: center

   Multilayer Signalling Model (MSM) overview. Segmented 3D cellular data across successive tissue layers (left) are used to construct a depth-resolved contact network. The MSM simulates lateral inhibition on this layered structure to predict SOP fate decisions (right), incorporating both apical and lateral cell-cell interactions. Adapted from Paci et al. (2026).

Citation
===========

If you use this model or documentation in your work, please cite:

- Paci, G., Berkemeier, F., Baum, B., Page, K. M., & Mao, Y. 3D epithelial cell topology tunes signaling range to promote precise patterning. *Proceedings of the National Academy of Sciences* 123(19): e2522727123 (2026). `https://www.pnas.org/doi/10.1073/pnas.2522727123 <https://www.pnas.org/doi/10.1073/pnas.2522727123>`_

Bugs, Questions and Comments
============================

If you encounter any issues or have questions about how to use the software, please open a `GitHub issue <https://github.com/fberkemeier/MultiLayer-NotchDelta/issues>`_. Feedback is a valuable part of the ongoing development. For comments or suggestions regarding improvements to the software or the documentation, you may also contact Francisco Berkemeier by email at fp409@cam.ac.uk.
