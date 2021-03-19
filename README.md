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


compilar:
./Tensor --> make
./FunctionSpace/OneDimensional/_OneD --> make
./FiniteElements/Assembly/_Assembly_ --> make cython_assembler_build
./FiniteElements/LocalAssembly/_KinematicsMeasures_ --> make
./VariationalPrinciple/_GeometricStiffness_ --> make
./VariationalPrinciple/_VolumetricStiffness_ --> make
./VariationalPrinciple/_ConstitutiveStiffness_ --> make
./MaterialLibrary/LLDispatch --> make MATERIAL=material_name
./FiniteElements/Assembly/_Assembly_ --> make ASSEMBLY_NAME=assembler_name

- Time integrator schema, explain.
