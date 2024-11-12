# Stochastic Modeling of Developmental Biology Concepts

This project is a MATLAB implementation of various stochastic models to simulate developmental biology concepts, with a focus on reaction-diffusion systems, compartmental modeling, and spatial heterogeneity in cellular processes. This was part of an undergraduate course project with professor David Green.

---

## Table of Contents
1. [Part A: Defining Concepts](#part-a-defining-concepts)
   - [1. Eukaryotic Cells, Maternal Effect Genes, and Somitogenesis](#1-eukaryotic-cells-maternal-effect-genes-and-somitogenesis)
   - [2. Sources of Spatial Heterogeneity in Bacterial Cells](#2-sources-of-spatial-heterogeneity-in-bacterial-cells)
   - [3. Skin Patterning in Zebrafish](#3-skin-patterning-in-zebrafish)
   - [4. Clock Wavefront Model](#4-clock-wavefront-model)
   - [5. Key Diffusion and Compartment Models](#5-key-diffusion-and-compartment-models)
   - [6. Reaction-Diffusion Model Using Finite Differences](#6-reaction-diffusion-model-using-finite-differences)
   - [7. Particle-Based Stochastic Reaction-Diffusion Model](#7-particle-based-stochastic-reaction-diffusion-model)
2. [Part B: Applying Concepts](#part-b-applying-concepts)
   - [1. Clock Oscillations and Somite Formation](#1-clock-oscillations-and-somite-formation)
   - [2. Two-Compartment Model for RNA and Protein Synthesis](#2-two-compartment-model-for-rna-and-protein-synthesis)
   - [3. Modeling Bacterial Growth in a Biofilm](#3-modeling-bacterial-growth-in-a-biofilm)
3. [Project Notes](#project-notes)

---

## Part A: Defining Concepts

### 1. Eukaryotic Cells, Maternal Effect Genes, and Somitogenesis
- **Eukaryotic Cells**: Multicellular organisms like plants, fungi, and animals are made up of eukaryotic cells, which have distinct compartments (organelles) with varying molecular concentrations necessary for biological processes.
- **Maternal Effect Genes**: These genes form gradients in the embryo, facilitating protein transcription for early developmental processes, like skin patterning.
- **Somitogenesis**: A segmentation process in vertebrate embryos where a molecular "clock" locks somites based on oscillating states, coordinating spinal cord formation.

### 2. Sources of Spatial Heterogeneity in Bacterial Cells
- **Examples**: Spatial heterogeneity in bacterial cells arises from localized concentrations of genetic material (nucleoids), ribosomes, membrane-bound vacuoles, and cell walls.

### 3. Skin Patterning in Zebrafish
- **Negative Feedback Loop**: Skin patterning in zebrafish is controlled by transcription factors in an activation-inhibition cycle, leading to spatial waves that define skin patterns based on differential diffusion rates.

### 4. Clock Wavefront Model
- **Activator-Inhibitor Pair**: Temporal oscillations are controlled by a pair of transcription factors, which propagate along the posterior to anterior wavefront, stabilizing states to define segmentation patterns.

### 5. Key Diffusion and Compartment Models
- **Fick’s Law of Diffusion**: Defines concentration change rates as proportional to the concentration gradient.
- **Compartment Model**: Models separate compartments in eukaryotic cells using Fick's First Law.
- **Brownian Motion**: Describes random particle movement in a medium, with equal probability in all directions.

### 6. Reaction-Diffusion Model Using Finite Differences
- **Simulation Process**: Finite differences and Forward Euler methods are used to simulate reaction-diffusion systems, approximating molecular concentrations across compartments and time points.

### 7. Particle-Based Stochastic Reaction-Diffusion Model
- **Model Structure**: This model uses molecule counts instead of concentrations, assigning reaction probabilities based on molecule availability. This yields stochastic variability, especially relevant for low-molecule environments.

## Part B: Applying Concepts

### 1. Clock Oscillations and Somite Formation
- **Wave Period and Growth**: The period of somite formation (X(d)) is influenced by clock oscillations and growth rate. Faster growth leads to larger somite sizes, as the wavefront locks the state of the clock in each segment.

### 2. Two-Compartment Model for RNA and Protein Synthesis
- **Flow Dynamics**: In this model, RNA synthesis occurs in the nucleus and protein synthesis in the cytoplasm. Concentration flows between compartments depend on concentration gradients and flow constants.

### 3. Modeling Bacterial Growth in a Biofilm
- **Quorum Sensing**: A compartment model where species X and Y represent signaling molecules, creating positive feedback for biofilm formation within bacteria. X induces aggregation in the biofilm, while Y suppresses it outside, establishing a spatial gradient.

## Project Notes
This project was completed independently, with references to course materials and online resources. MATLAB was used for the simulations, based on lecture slides and background literature.

---
