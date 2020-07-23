# vdp

This is a C++ implementation of the distortion point detection method in the following paper:

[Shuangming Chai](https://kfckfckf.github.io/), [Xiao-Ming Fu](http://staff.ustc.edu.cn/~fuxm), and [Ligang Liu](http://staff.ustc.edu.cn/~lgliu).
Voting for Distortion Points in Geometric Processing.
IEEE Transactions on Visualization and Computer Graphics, 2019, Early Access.
DOI: [10.1109/TVCG.2019.2947420](https://doi.org/10.1109/TVCG.2019.2947420)

The code is written by Shuangming Chai using Microsoft Visual Studio 2019.

## External Libraries

All the required libraries are the latest versions.

* [Eigen](http://eigen.tuxfamily.org/)
* [OpenMesh](https://www.openmesh.org/)
* [PARDISO](https://www.pardiso-project.org/) (academic license) or [Intel MKL](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html)
* OpenMP

## Usage

```
vdp <mesh> [-a] [-g] [-p <int>] [-q <int>] [-r <int>] [-s <int>] [-t <int>] [-v]
```

Command arguments:

* `mesh`: a genus zero mesh (.obj .off .stl .ply).
* `-a`/`--intermediate`: save intermediate results. (optional)
* `-g`/`--geodesic`: use geodesic distance. (optional)
* `-p`/`--nsimplify`: vertex number of simplified mesh. (optional) default: 13000
* `-q`/`--nrandom`: number of random cutting process. (optional) default: 10
* `-r`/`--nvoting`: voting threshold. (optional) default: 2
* `-s`/`--nlocal`: specify the feature size N. (optional) default: 13
* `-t`/`--npost`: n-ring neighborhood. (optional) default: 5
* `-v`/`--verbose`: output debug info. (optional)

### Input
A closed genus-zero surface mesh.

### Output
Suppose that the input file name is XXX.obj, there are three output files:

* XXX_cut.txt: text file including vertex indices and edge indices.
* XXX_landmarks.txt: text file including vertex indices of the distortion points.
* XXX_open.obj: the open mesh cut by our method.

Note: The indices of vertices and edges are on the original mesh, and they can be changed on the resulting open mesh.

## Example

Example 1:
```
vdp bear.obj -a -s 13
```

Example 2:
```
vdp test.obj --verbose --nsimplify 10000 --nlocal 10
```

## Changelog

### Version 0.2 (2020-07-23)
* Add Intel MKL Pardiso solver support.
* Add a class ParameterManager to manage parameters.
* Add verbose parameter to control output debug messages.
* Reformat the code with a consistent coding style.
* Other improvements and fixed bugs.

### Version 0.1 (2019-10-17)
* Initial version
