# Computational Neuroscience — EC60007

Original MATLAB implementations and the report for the Computational Neuroscience course (EC60007), IIT Kharagpur.

## Table of Contents

- [Overview](#overview)
- [Contents](#contents)
- [Concepts Covered](#concepts-covered)
- [Tools](#tools)
- [Repository Structure](#repository-structure)
- [Running the Code](#running-the-code)
- [Academic Context](#academic-context)

## Overview

| Item | Detail |
|------|--------|
| Course | EC60007 — Computational Neuroscience |
| Author | Jai Raj (23EC10031) |
| Institution | IIT Kharagpur, Department of E&ECE |
| Instructor | Sharba Bandyopadhyay |
| Work | Project II: Morris-Lecar and Hodgkin-Huxley neuron models |

## Contents

| Project | Topic | Files |
|---------|-------|-------|
| [Morris-Lecar and Hodgkin-Huxley](computational-neuroscience/morris-lecar-hodgkin-huxley/) | Equilibrium, phase plane, nullclines, threshold, firing rate, bistability, anode break | 2 MATLAB scripts, 27 figures, report |

## Concepts Covered

- Morris-Lecar model: equilibrium, vector field, nullclines, phase-plane analysis.
- Threshold behaviour and firing-rate curves.
- Bistability and anode-break excitation.
- Hodgkin-Huxley model: gating variables, threshold scans, full vs reduced dynamics.
- Numerical integration with `ode45` and tight tolerances.

## Tools

- MATLAB (base, no toolboxes beyond `ode45`/`fzero`).
- Report as PDF.

## Repository Structure

```text
computational-neuroscience/
├── README.md
├── .gitignore
└── computational-neuroscience/
    ├── README.md
    └── morris-lecar-hodgkin-huxley/
        ├── README.md
        ├── report.pdf
        ├── src/
        │   ├── project02.m
        │   └── project_ii.m
        ├── figures/ (27 PNGs)
        └── supplementary/
            ├── project_i_instructions.pdf
            └── project_ii_instructions.pdf
```

## Running the Code

```matlab
cd morris-lecar-hodgkin-huxley/src
project_ii   % refined version (recommended)
% or
project02    % earlier version
```

Both scripts generate the figures saved in `figures/`. `project_ii.m` is the refined version; its header notes that AI was used only for formatting/refinement, not for the underlying analysis.

## Academic Context

Individual coursework by Jai Raj (23EC10031). An empty `project01.m` and the NPTEL transcripts / past-year questions from the same course are intentionally not included. No license is applied.
