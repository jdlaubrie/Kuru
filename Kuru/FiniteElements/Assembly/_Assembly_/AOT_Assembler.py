from __future__ import print_function
import os, platform, sys, subprocess

__all__ = ['AOTConfigure','AOTClean']



def MaterialList():

    ll_materials_mech          = [  "_NeoHookean_2_",
                                    "_AnisotropicFungQuadratic_",
                                    "_ArterialWallMixture_",
                                    ]

    return ll_materials_mech

def execute(_cmd):
    _process = subprocess.Popen(_cmd, shell=True)
    _process.wait()



def AOTConfigure():

    ll_materials_mech = MaterialList()


    for material in ll_materials_mech:
        material_specific_assembler = "_LowLevelAssemblyDF_" + material

        header_f = material_specific_assembler + ".h"
        cython_f = material_specific_assembler + ".pyx"

        execute("cp _LowLevelAssemblyDF_.h " + header_f)
        execute("cp _LowLevelAssemblyDF_.pyx " + cython_f)

        # Read the header file
        f = open(header_f, "r")
        contents_h = f.readlines()
        f.close()

        # Modify header file
        contents_h[0] = "#ifndef " + material_specific_assembler.upper() + "_H\n"
        contents_h[1] = "#define " + material_specific_assembler.upper() + "_H\n"
        contents_h[5] = '#include "' + material + '.h"\n'

        for counter, line in enumerate(contents_h):
            if "_GlobalAssemblyDF_" in line:
                contents_h[counter] = contents_h[counter].replace("_GlobalAssemblyDF_", "_GlobalAssemblyDF_"+material)

            if "auto mat_obj" in line:
                if material == "_NeoHookean_" or material == "_LinearElastic_":
                    contents_h[counter] = "    auto mat_obj = " + material + "<Real>(mu,lamb);\n"
                elif material == "_MooneyRivlin_":
                    contents_h[counter] = "    auto mat_obj = " + material + "<Real>(mu1,mu2,lamb);\n"
                elif material == "_NearlyIncompressibleMooneyRivlin_":
                    contents_h[counter] = "    auto mat_obj = " + material + "<Real>(mu1,mu2,lamb);\n"
                elif material == "_AnisotropicMooneyRivlin_1_":
                    contents_h[counter] = "    auto mat_obj = " + material + "<Real>(mu1,mu2,mu3,lamb);\n"

        if material == "_AnisotropicMooneyRivlin_1_":
            for counter, line in enumerate(contents_h):
                if "mat_obj.KineticMeasures" in line:
                    contents_h[counter] = "        mat_obj.KineticMeasures(stress, hessian, ndim, ngauss, F, &anisotropic_orientations[elem*ndim]);\n"

        # Turn off geometric stiffness
        if material == "_LinearElastic_":
            for counter, line in enumerate(contents_h):
                if "std::fill" in line and "geometric_stiffness" in line:
                    contents_h[counter] = ""
                # DANGEROUS - this is relative to current line - we delete 13 lines down
                if "_GeometricStiffnessFiller_" in line:
                    contents_h[counter:counter+13] = ""


        contents_h[-1] = "#endif // " + material_specific_assembler.upper() + "_H"

        # Write
        f = open(header_f, 'w')
        for item in contents_h:
            f.write(item)


        # Read the cython wrapper file
        f = open(cython_f, "r")
        contents_c = f.readlines()
        f.close()


        for counter, line in enumerate(contents_c):
            if "_LowLevelAssemblyDF_" in line:
                contents_c[counter] = contents_c[counter].replace("_LowLevelAssemblyDF_", "_LowLevelAssemblyDF_" + material)
            if "_GlobalAssemblyDF_" in line:
                contents_c[counter] = contents_c[counter].replace("_GlobalAssemblyDF_", "_GlobalAssemblyDF_" + material)

            if "mu1, mu2, lamb =" in line:
                if material == "_NeoHookean_" or material == "_LinearElastic_":
                    contents_c[counter] = "    mu, lamb = material.mu, material.lamb\n"
                elif material == "_MooneyRivlin_":
                    contents_c[counter] = "    mu1, mu2, lamb = material.mu1, material.mu2, material.lamb\n"
                elif material == "_NearlyIncompressibleMooneyRivlin_":
                    contents_c[counter] = "    mu1, mu2, lamb = material.mu1, material.mu2, material.lamb\n"
                elif material == "_AnisotropicMooneyRivlin_1_":
                    contents_c[counter] = "    mu1, mu2, mu3, lamb = material.mu1, material.mu2, material.mu3, material.lamb\n"

        # Write
        f = open(cython_f, 'w')
        for item in contents_c:
            f.write(item)


def AOTClean():

    ll_materials_mech = MaterialList()

    for material in ll_materials_mech:
        material_specific_assembler = "_LLADF_" + material

        header_f = material_specific_assembler + ".h"
        cython_f = material_specific_assembler + ".pyx"

        execute("rm -rf " + header_f)
        execute("rm -rf " + cython_f)


if __name__ == "__main__":

    args = sys.argv

    _op = None

    if len(args) > 1:
        for arg in args:
            if arg == "clean" or arg == "configure":
                if _op is not None:
                    raise RuntimeError("Multiple conflicting arguments passed to setup")
                _op = arg

    if _op == "clean":
        AOTClean()
    elif _op == "configure":
        AOTConfigure()
    else:
        AOTClean()
        AOTConfigure()

