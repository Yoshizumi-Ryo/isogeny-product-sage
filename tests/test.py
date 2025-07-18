

######################################################################
# In this module, you can compute examples of (ell,ell)-isogeny using.
######################################################################

import pytest
from sage.all import GF,is_prime,gcd,EllipticCurve
from src.func_chain import Total_computation
from src.func_E0 import Pre_Random_Isog_Images
from src.func_test import Average_time

#base field======================================================

# PYTHONPATH=src sage -python -m pytest -s tests/test.py



@pytest.fixture
def setup_isogeny():
    p = 276154505650672190920223
    D = 3**6 * 7 * 11**4 * 19 * 29**3 * 67 * 79
    d = 4632525325
    Fp2=GF(p**2)
    assert(is_prime(p))
    assert(p%4==3)
    assert((p+1)%D==0)
    assert(D>d)
    assert(gcd(D,d)==1)
    assert((D-d)*d>p)
    E0 = EllipticCurve(Fp2, [0, 0, 0, +1, 0])#E0:y^2=x^3+x
    assert(E0.is_supersingular())
    assert(E0.order()==(p+1)**2)
    assert(E0.j_invariant()==1728)
    E1=E0
    E2=E0
    [P1,P2],[Q1,Q2]=Pre_Random_Isog_Images(p,D,d,E0)
    assert P1.weil_pairing(Q1,D)*P2.weil_pairing(Q2,D)==1 
    #E_1
    ord1=19
    x1_E1,_=E1.torsion_basis(ord1)
    ord2=7
    x2_E1,_=E1.torsion_basis(ord2)
    return p,D,d,Fp2,E1,E2,P1,P2,Q1,Q2,ord1,x1_E1,ord2,x2_E1




def test_computation(setup_isogeny):
    p,D,d,Fp2,E1,E2,P1,P2,Q1,Q2,ord1,x1_E1,ord2,x2_E1=setup_isogeny
    (E1_cod,E2_cod,pt_fx_E1_list,pt_fx_E2_list),_,_=Total_computation(
    E1,E2,P1,P2,Q1,Q2,[x1_E1,x2_E1],D,"One","Existing")
    assert pt_fx_E1_list[0].order()==ord1
    assert pt_fx_E1_list[1].order()==ord2
    assert pt_fx_E2_list[0].order()==ord1
    assert pt_fx_E2_list[1].order()==ord2



# Average_time(10,1,3)
# Average_time(10,2,3)
