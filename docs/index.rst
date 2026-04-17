.. _index:

##########################
StACKER documentation
##########################

StACKER manipulates the outputs of an MD simulation and analyzes 
the pi-stacking interactions. Creates a "Pi-Stacking Fingerprint" 
for a structure at each frame. Presents Pi-stacking interactions between 
two residues through analysis of their relevant movement.

.. grid:: 1 1 2 2
    :gutter: 2 3 4 4

    .. grid-item-card::
        :text-align: center

        .. raw:: html

            <img src="_static/index_images/getting_started.svg" style="height: 60px; width: 60px; margin-left: auto; margin-right: auto; margin-top: 10px;" />

        Getting Started
        ^^^

        Learn how to install StACKER and
        prepare Molecular Dynamics Data

        +++

        .. button-ref:: getting_started
            :expand:
            :color: secondary
            :click-parent:

            To the installation instructions

    .. grid-item-card::
        :text-align: center

        .. raw:: html

            <img src="_static/index_images/python_logo.webp" style="height: 60px; width: 60px; margin-left: auto; margin-right: auto; margin-top: 10px;" />

        Python Modules
        ^^^

        See detailed instructions on StACKER's 
        functions and modules when using Python

        +++

        .. button-ref:: modules
            :expand:
            :color: secondary
            :click-parent:

            To the Python user guide

    .. grid-item-card::
        :text-align: center

        .. raw:: html

            <img src="_static/index_images/system_stacking_fingerprint.png" style="height: 200px; width: 200px; margin-left: auto; margin-right: auto; margin-top: 10px;" />

        Create System Stacking Fingerprints (SSFs)
        ^^^

        Pipeline to create System Stacking Fingerprints (SSFs)
        using the Python Options

        +++

        .. button-ref:: system_stacking
            :expand:
            :color: secondary
            :click-parent:

            To the SSF creation guide

    .. grid-item-card::
        :text-align: center

        .. raw:: html

            <img src="_static/index_images/pairwise_stacking_fingerprint.png" style="height: 200px; width: 250px; margin-left: auto; margin-right: auto; margin-top: 10px;" />

        Create Pairwise Stacking Fingerprints (PSFs)
        ^^^

        Pipeline to create Pairwise Stacking Fingerprints (PSFs)
        using the Python Options

        +++

        .. button-ref:: pairwise_stacking
            :expand:
            :color: secondary
            :click-parent:

            To the PSF creation guide

    .. grid-item-card::
        :text-align: center

        .. raw:: html

            <img src="_static/index_images/command_line.png" style="height: 60px; width: 60px; margin-left: auto; margin-right: auto; margin-top: 10px;" />

        Command Line Options
        ^^^

        See a detailed pipeline to generate SSFs
        and PSFs using the Command Line Options
        +++

        .. button-ref:: command_line_options
            :expand:
            :color: secondary
            :click-parent:

            To the Command Line user guide

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   Getting Started <getting_started>
   API Documentation <modules>
   Create System Stacking Fingerprints (SSFs) <system_stacking>
   Create Pairwise Stacking Fingerprints (PSFs) <pairwise_stacking>
   Command Line Options <command_line_options>
   Molecular Dynamics <molecular_dynamics>

