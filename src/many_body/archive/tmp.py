from triqs.gf.meshes import MeshDLRImFreq
beta, wmax, eps = 100.0, 2.0, 1e-8
m = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=wmax, eps=eps)
print(len(m))

