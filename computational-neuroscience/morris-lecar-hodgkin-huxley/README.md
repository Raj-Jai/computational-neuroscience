# Morris-Lecar and Hodgkin-Huxley Neuron Models

**Course:** EC60007 Computational Neuroscience · **Author:** Jai Raj (23EC10031) · **Instructor:** Sharba Bandyopadhyay

[Read the report (PDF)](report.pdf)

## Objective

Analyse single-neuron dynamics for the Morris-Lecar (ML) and Hodgkin-Huxley (HH) models: equilibrium and phase-plane structure, threshold behaviour, firing-rate curves, bistability and anode-break excitation.

## Implementation

- `project02.m` — earlier version of the analysis.
- `project_ii.m` — refined version (recommended); uses `ode45` with tight tolerances and `fzero` for equilibria. Header notes that AI was used only for formatting/refinement.

Both scripts generate the figures in `figures/` (Q02–Q18).

## Files

| File | Description |
|------|-------------|
| `report.pdf` | Full report (Project II) |
| `src/project02.m` | Earlier MATLAB implementation |
| `src/project_ii.m` | Refined MATLAB implementation |
| `figures/` | 27 plots (phase planes, nullclines, firing rates, threshold scans, HH gating variables) |
| `supplementary/` | Project I and II instruction sheets |

## Running

```matlab
cd src
project_ii
```

Requires base MATLAB only.
