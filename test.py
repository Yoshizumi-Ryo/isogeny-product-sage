

######################################################################
# In this module, you can compute examples of (ell,ell)-isogeny using.
######################################################################


from sage.all import GF,is_prime,gcd,EllipticCurve
from func_chain import Preparation_dim1_theta,Odd_deg_isogeny_chain,Splitting,Total_computation
from func_E0 import Pre_Random_Isog_Images
#base field======================================================


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


#elliptic curve===================================================


#E0:y^2=x^3+x
E0 = EllipticCurve(Fp2, [0, 0, 0, +1, 0])

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




#compute theta coordinate in dim 1 =========================

tnp_E1,first_basis_E1,ext_P1,ext_Q1,ext_x_E1_list,tnp_E2,first_basis_E2,ext_P2,ext_Q2=Preparation_dim1_theta(E1,E2,P1,P2,Q1,Q2,[x1_E1],D)



tnp_cod,tc_x_cod_list,_=Odd_deg_isogeny_chain(
    tnp_E1,first_basis_E1,ext_P1,ext_Q1,ext_x_E1_list,tnp_E2,first_basis_E2,ext_P2,ext_Q2,D,"Square","Proposed")



E1_cod,E2_cod,pt_fx_E1_list,pt_fx_E2_list=Splitting(tnp_cod,tc_x_cod_list)



#==================================================


E1_cod,E2_cod,pt_fx_E1_list,pt_fx_E2_list=Total_computation(
    E1,E2,P1,P2,Q1,Q2,[x1_E1,x2_E1],D,"One","Existing")


assert pt_fx_E1_list[0].order()==ord1
assert pt_fx_E1_list[1].order()==ord2
assert pt_fx_E2_list[0].order()==ord1
assert pt_fx_E2_list[1].order()==ord2



#==================================================


from func_test import Average_time

Average_time(10,1,3)
Average_time(10,2,3)

