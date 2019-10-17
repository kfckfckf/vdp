vdp
===

This is a C++ implementation of the distortion point detection method in "Voting Distortion Points for Geometric Processing".
The code is compiled using Microsoft Visual Studio 2017 with OpenMP support for parallel programming.

Library
---

* [Eigen](http://eigen.tuxfamily.org/)
* [OpenMesh](https://www.openmesh.org/)
* [PARDISO](https://www.pardiso-project.org/)
* OpenMP support

Usage
---

```
vdp <mesh> [-a] [-s N]
```

Command arguments:
* mesh: a genus zero mesh.
* -a: save intermediate results. (optional)
* -s: specify the feature size N. (optional)

Example
---

```
vdp bear.obj -a -s 13
```
