# StACKER 2.0
**St**acking **A**nalysis and **C**onformational **K**inetics for **E**xamining **R**esidues

Developed by Pete Hwang ([phwang@wesleyan.edu](mailto:phwang@wesleyan.edu)) in the [Weir Lab](https://weirlab.research.wesleyan.edu/) at Wesleyan University.

StACKER 2.0 is a Python toolkit for systems-level π-stacking analysis of molecular dynamics (MD) trajectories. It generates **System Stacking Fingerprints (SSFs)** for global residue–residue stacking geometry and **Pairwise Stacking Fingerprints (PSFs)** for detailed analysis of specific residue pairs.

Originally developed by Eric Sakkas in the Weir Lab at Wesleyan University, StACKER has been expanded in StACKER 2.0 to support **mixed RNA–protein systems**, including **nucleotides, aromatic amino acids, and arginine**.

## Project Overview

StACKER analyzes π-stacking interactions across biomolecular structures and trajectories by representing residue–residue geometry in an interpretable fingerprint framework.

The package is centered around two complementary representations:

- **System Stacking Fingerprints (SSFs)**  
  Symmetric residue-by-residue distance matrices that summarize center-of-geometry relationships across an entire structure or trajectory.

- **Pairwise Stacking Fingerprints (PSFs)**  
  Bottaro-style polar-coordinate plots that describe how one residue moves relative to another across a trajectory.

Together, these tools allow users to study stacking behavior at both the **systems level** and the **individual residue-pair level**.

## What's New in StACKER 2.0

Compared with the original StACKER framework, StACKER 2.0 introduces:

- support for **aromatic amino acids and arginine** in addition to nucleotides
- updated default center-of-geometry definitions for:
  - pyrimidines
  - purines
  - phenylalanine
  - tyrosine
  - tryptophan
  - histidine
  - arginine guanidinium groups
- amino-acid-compatible **Pairwise Stacking Fingerprints (PSFs)**
- class-aware default cutoffs for comparison workflows
- substantial SSF performance improvements through a redesigned and parallelized computation pipeline

## Supported Residue Classes

StACKER 2.0 currently supports stacking analysis for:

### Nucleotides
- A
- G
- C
- U
- T
- inosine
- common terminal variants included in trajectory/topology naming conventions

### Aromatic amino acids
- Phe
- Tyr
- Trp
- His  
  including common protonation-state variants:
  - HID
  - HIE
  - HIP

### Charged planar group
- Arg  
  using the guanidinium heavy atoms as the default stacking-relevant group

## Installation

Clone the repository and install StACKER 2.0 locally:

```bash
git clone https://github.com/phwang25/stacker2.git
cd stacker2
pip install -e .