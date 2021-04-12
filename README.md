# kuru

This is a Finite-Element code to solve my tesis problem about Growth and Remodeling of aneurysms on the Aorta artery. The code is mainly base on Florence.

- New follower load (pressure) in Assembly module.
- New ArterialWallMixture in Material Library.
- New Explicit Growth and Remodeling time integrator.

Prerequisites:

- Fastor
- Cython
- NumPy
- SciPy
- OpenBLAS


To compile the project:

```
$ python setup.py build
```

To clean sources and compilation files:

```
$ python setup.py clean
```
