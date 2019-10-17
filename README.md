vdp
===

This is a C++ implementation of the distortion point detection method in the following paper:
*Voting Distortion Points for Geometric Processing*
[Shuangming Chai](https://kfckfckf.github.io/), [Xiao-Ming Fu](http://staff.ustc.edu.cn/~fuxm), and [Ligang Liu](http://staff.ustc.edu.cn/~lgliu).
IEEE Transactions on Visualization and Computer Graphics, 2019.
The code is written by Shuangming Chai using Microsoft Visual Studio 2017 with OpenMP.

External Libraries
---

All the required libraries are the latest versions.
* [Eigen](http://eigen.tuxfamily.org/)
* [OpenMesh](https://www.openmesh.org/)
* [PARDISO](https://www.pardiso-project.org/) (academic license)
* OpenMP

Usage
---

```
vdp <mesh> [-a] [-g] [-p <int>] [-q <int>] [-r <int>] [-s <int>] [-t <int>]
```

Command arguments:
* mesh: a genus zero mesh.
* -a: save intermediate results. (optional)
* -g: use geodesic distance. (optional)
* -p: vertex number of simplified mesh. (optional) default: 13000
* -q: number of random cutting process. (optional) default: 10
* -r: voting threshold. (optional) default: 2
* -s: specify the feature size N. (optional) default: 13
* -t: n-ring neighborhood. (optional) default: 5

Example
---

```
vdp bear.obj -a -s 13
```
