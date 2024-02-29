from pyscf import gto, scf
from scipy.linalg import fractional_matrix_power
import numpy as np

"""mol = gto.M(atom = F'O 0.0 0.0 0.0; H 1.0 0.0 0.0', basis = 'ccpvdz', charge = -1)
mol.nelec = (5, 4)
mf = scf.RHF(mol)
uhf = scf.UHF(mol)
E = uhf.kernel()
SCF_E_pyscf = mf.kernel()
print(mf.get_occ())
print(mf.analyze())
print(mf.get_hcore())
print (mf.get_init_guess())
e_neg = -75.3297749390419
e_neu_1  =-75.3884393376133"""
"""mol = gto.M(atom = F'H 1.0 0.0 0.0', basis = 'ccpvdz', charge = 1)
mol.nelec = (1, 0)
mf = scf.RHF(mol)
uhf = scf.UHF(mol)
E = uhf.kernel()
SCF_E_pyscf = mf.kernel()
e_neu_2 = -0.499278403419583
e_pos = 0
print(mf.get_occ())
print(mf.analyze())
print(mf.get_hcore())
print(e_neg + e_pos)
print (e_neu_1 + e_neu_2)"""
mol = gto.M(atom = F'O 0.0 0.0 0.0; H 1.0 0.0 0.0 ; H 0.0 300.0 0.0', basis = 'ccpvdz')
mf = scf.RHF(mol)
uhf = scf.UHF(mol)
E = uhf.kernel()
SCF_E_pyscf = mf.kernel()
print(mf.get_occ())
print(mf.analyze())

