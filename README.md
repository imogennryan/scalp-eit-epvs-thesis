# Scalp EIT head model for enlarged perivascular spaces

Code for my Biomedical Engineering Honours thesis at the University of Sydney (2026):

**Computational Feasibility of Scalp Electrical Impedance Tomography for Monitoring Glymphatic Dysfunction via Detection of Enlarged Perivascular Spaces**

Imogen Ryan. Supervisors: Prof. Alistair McEwan and Monzurul Alam.

## What it does

The code simulates whether a build-up of fluid around the small blood vessels of the brain (enlarged perivascular spaces) changes the voltages measured on the scalp by more than the noise of an EIT machine. Everything is simulated. There is no patient data.

1. **Aim 1. Build and validate the head model.** A four-layer head model (scalp, skull, CSF, brain) is built in EIDORS from the Tizzard adult head mesh, with electrodes at the 10-20 and 10-10 EEG positions. Eight tests check it.
2. **Aim 2. Simulate the readings.** A region of extra fluid is added to the brain at different scalp points, depths, fluid amounts, shapes and electrode layouts (1,578 simulations). The change in the scalp voltages is compared with instrument noise.

This is a forward model only. It asks whether a difference can be detected, not where it is.

## What you need

- MATLAB (written and tested on R2026a)
- EIDORS 3.12, from https://eidors3d.sourceforge.net
- Python 3 with `openpyxl` (only used to colour one spreadsheet)

## How to run it

In MATLAB:

```matlab
run('<your EIDORS folder>/eidors/startup.m')   % start EIDORS
run('run_all.m')                                % run everything in order
```

This takes about 13 minutes on a laptop. It writes everything into this folder:

| Folder | What is in it |
|---|---|
| `models/` | the electrode layouts and the healthy head for each layout |
| `data/` | the change in the readings for every simulation |
| `results/` | the tables of results (CSV) |
| `figures/` | the figures in the thesis (Figures 3.1 to 3.8 and 4.1 to 4.7) |

Each script can also be run on its own once the scripts before it in `run_all.m` have run.

## Folders

| Folder | What is in it |
|---|---|
| `code/aim1_head_model/` | building and testing the head model |
| `code/aim2_simulations/` | the fluid simulations, noise and Results |
| `code/tools/` | a small Python script that colours one spreadsheet |
| `models/tizzard_meshes/` | the seven Tizzard adult head meshes |

Each file is described in Appendix A of the thesis (Table A.1), and each script starts with a short description of what it does.

## Head meshes

The head meshes are by Andrew Tizzard and are distributed with EIDORS under a Creative Commons licence that requires attribution. If you use them, please cite:

> Tizzard, A. and Bayford, R.H. (2007). Improving the Finite Element Forward Model of the Human Head by Warping using Elastic Deformation. Physiol. Meas. 28:S163-S182.

See `models/tizzard_meshes/README.md`.

## Use of AI

I wrote this code with the help of Anthropic's Claude (Claude Code), as declared in the Generative AI Use Statement of the thesis. I ran every script and checked every result.

## Licence

The code is released under the MIT licence (see `LICENSE`). This does not cover the Tizzard meshes or EIDORS, which have their own licences.
