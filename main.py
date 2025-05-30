

##################################################################
# In this module, you can compute random (ell,ell)-isogeny from  E0*E0 for fixed p,l.
# p: characteristic
# l: degree
# E0: y^2=x^3+x over F_{p^2}.
# For more details, see "test.py" in the same package.
# 
# For example, write $ sage main.py "Proposed" "One"  1
################################################################



if __name__ == "__main__":
    import sys
    from sage.all import GF,is_prime,gcd,EllipticCurve
    from func_chain import Total_computation
    from func_E0 import Pre_Random_Isog_Images
    args = sys.argv
    assert(len(args) == 4)
    _,Proposed_or_Exisited,One_or_Square,parameter= args
    assert(Proposed_or_Exisited in {"Proposed","Existing"})
    assert(One_or_Square in {"One","Square"})
    assert(parameter in {"1","2"})
    #base field======================================================
    if parameter=="1":
        p = 276154505650672190920223
        D = 3**6 * 7 * 11**4 * 19 * 29**3 * 67 * 79
        d = 2**(30)
        ord_1,ord_2 = 7,19
    else: # parameter=="2"
        p=0x1935BECE108DC6C0AAD0712181BB1A414E6A8AAA6B510FC29826190FE7EDA80F
        D=3*7**(16)*17**9*31**8*311*571*1321
        d=2**(150)
        ord_1,ord_2 =17,31
    print("p=",p)
    Fp2=GF(p**2)
    assert(is_prime(p))
    assert(p%4==3)
    assert((p+1)%D==0)
    assert(D>d)
    assert(gcd(D,d)==1)
    assert((D-d)*d>p)
    #elliptic curve===========================
    #E0:y^2=x^3+x
    E0 = EllipticCurve(Fp2, [0, 0, 0, +1, 0])
    assert(E0.is_supersingular())
    assert(E0.order()==(p+1)**2)
    assert(E0.j_invariant()==1728)
    E1=E0
    E2=E0
    [P1,P2],[Q1,Q2]=Pre_Random_Isog_Images(p,D,d,E0)
    assert P1.weil_pairing(Q1,D)*P2.weil_pairing(Q2,D)==1 
    #on E_1
    x1_E1,_=E1.torsion_basis(ord_1)
    x2_E1,_=E1.torsion_basis(ord_2)
    #computation ==========================
    E1_cod,E2_cod,pt_fx_E1_list,pt_fx_E2_list=Total_computation(
    E1,E2,P1,P2,Q1,Q2,[x1_E1,x2_E1],D,One_or_Square,Proposed_or_Exisited)



