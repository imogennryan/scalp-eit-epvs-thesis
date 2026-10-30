# Tizzard adult head meshes

The seven four-layer adult head meshes (scalp, skull, CSF, brain) used in this thesis. SAH262 is the mesh used for every result. The others are used for the mesh tests.

## Attribution

These meshes are distributed under a Creative Commons Artistic License with attribution. Any use must cite:

> Tizzard, A. and Bayford, R.H. (2007). *Improving the Finite Element Forward Model of the Human Head by Warping using Elastic Deformation.* Physiol. Meas. **28**:S163-S182.

The geometry is described in:

> Tizzard, A., Horesh, L., Yerworth, R.J., Holder, D.S. and Bayford, R.H. (2005). *Generating accurate finite element meshes for the forward model of the human head in EIT.* Physiol. Meas. **26**(2):S251-S261.

Author: Andrew Tizzard, Middlesex University. Source: https://eidors3d.sourceforge.net/data_contrib/at-head-mesh/

The MIT licence in this repository does not cover these meshes.

## Contents

| File | Elements | Nodes | SHA-256 (first 16) |
|---|---|---|---|
| SAH028.mat | 27 684 | 5 226 | 2b67e86da35d2c64 |
| SAH031.mat | 31 111 | 5 939 | 1544aa50d7506073 |
| SAH035.mat | 35 073 | 6 556 | 19f633533ec5f964 |
| SAH049.mat | 49 172 | 9 142 | fec903307729b3bf |
| SAH077.mat | 77 196 | 14 188 | c782ca96f5868821 |
| SAH132.mat | 131 672 | 24 688 | a00c64a9ce83cec8 |
| SAH262.mat | 262 368 | 48 394 | f2f70f70720f0123 |

Element and node counts match the published table on the EIDORS page (checked by `a1_verify_mesh_family.m`).

All seven are the same head at different mesh densities.

## File structure

Each `.mat` contains:

- `vtx` - nodal positions, in millimetres
- `tri` - tetrahedra
- `materials` - struct with `names`, `sigma`, `ref` (`ref` indexes each element
  into `names`)

## Tissue order is not the same in every file

| Mesh | Order in `materials.names` |
|---|---|
| SAH031 | **scalp \| skull \| CSF \| brain** |
| all others | skull \| brain \| scalp \| CSF |

So `ref == 2` is brain in six of the files but skull in SAH031. The code always finds a tissue by its name:

```matlab
names  = strtrim(cellstr(S.materials.names));
iBrain = find(strcmpi(names,'brain'));
```

## Units

Nodes are in millimetres and EIDORS works in metres. `a1_build_head_model.m` converts them when the mesh is loaded.

## Conductivities in the files

Each mesh file includes a set of conductivities with no documented source. This thesis does not use them. It uses scalp 0.23, skull 0.05, CSF 1.79 and brain 0.10 S/m (McDermott et al. 2020; Baumann et al. 1997), set in `a1_build_head_model.m`.
