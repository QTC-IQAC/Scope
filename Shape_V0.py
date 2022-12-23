#
#cellpath = "/Users/ycho/Projects/cell2mol_materialcloud/T-TMCs/1-Iron/CEMWUO_TMC_1.gmol"
#metal_dir = "1-Iron"
#with open(cellpath, "rb") as pickle_file:
#    mol = pickle.load(pickle_file)
#if mol.hapticity == False:
#    print(f"{mol.refcode} CN {mol.totmconnec}") # save TMC with no hapticity of ligands
#    filename = "/Users/ycho/Projects/cell2mol_materialcloud/T-TMCs/{}/coord_sphere/{}_shape.xyz".format(metal_dir, mol.refcode)
#    with open(filename, "w") as fil:
#        print(mol.totmconnec+1, file=fil)
#        print(mol.refcode, file=fil)
#        met = mol.metalist[0]
#        print("%s  %.6f  %.6f  %.6f" % (met.label, met.coord[0], met.coord[1], met.coord[2]),file=fil)
#        for lig in mol.ligandlist:
#            for cn, ml, mlc in zip(lig.mconnec, lig.labels, lig.coord):
#                if cn>= 1:
#                    print("%s  %.6f  %.6f  %.6f" % (ml, mlc[0], mlc[1], mlc[2]),file=fil)
#else:
#    pass
