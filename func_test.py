
from sage.all import GF,is_prime,gcd,EllipticCurve
from func_chain import Total_computation
from func_E0 import Pre_Random_Isog_Images


def Test(Proposed_or_Existing,One_or_Square,evaluation_pts=1,parameter=3):
    assert(Proposed_or_Existing in {"Proposed","Existing"})
    assert(One_or_Square in {"One","Square"})
    assert(evaluation_pts in {1,2})
    assert(parameter in {1,2,3})
    print(" ")
    print(Proposed_or_Existing,One_or_Square)
    #base field======================================================
    if parameter==1:
        p = 276154505650672190920223
        D = 3**6 * 7 * 11**4 * 19 * 29**3 * 67 * 79
        d = 2**(30)
        ord_1,ord_2=7,11
    elif parameter==2:
        p=0x1935BECE108DC6C0AAD0712181BB1A414E6A8AAA6B510FC29826190FE7EDA80F
        D=3*7**(16)*17**9*31**8*311*571*1321
        d=2**(150)
        ord_1,ord_2=17,31
    else:
        p = 110564446907225951023
        D = 7**3 * 11 * 19 * 23**2 * 53 * 137
        d=2**(37)
        ord_1,ord_2=7,11
    print("the_number_of_evaluation_points:",evaluation_pts)
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
    if evaluation_pts==1:
        pt_list=[x1_E1]
    else:
        pt_list=[x1_E1,x2_E1]
    #computation ==========================
    _,time1,time2=Total_computation(
        E1,E2,P1,P2,Q1,Q2,pt_list,D,One_or_Square,Proposed_or_Existing)
    print(" ")
    return time1,time2



def All_tests(evaluation_pts=1,parameter=3):
    time1_ps,time2_ps=Test("Proposed","Square",evaluation_pts,parameter)
    time1_po,time2_po=Test("Proposed","One"   ,evaluation_pts,parameter)
    time1_es,time2_es=Test("Existing","Square",evaluation_pts,parameter)
    time1_eo,time2_eo=Test("Existing","One"   ,evaluation_pts,parameter)
    return [time1_ps,time2_ps,
           time1_po,time2_po,
           time1_es,time2_es,
           time1_eo,time2_eo]



def Average_time(n:int,evaluation_pts=1,parameter=3):
    """  
    n: the numerator of aveage.
    evaluation_pts: the number of evaluation points of the form (x1,0).
    parameter: 1,2,3.
    """
    sum_time=[0,0,0,0,0,0,0,0]
    for i in range(n):
        this_time=All_tests(evaluation_pts,parameter)
        for t in range(8):
            sum_time[t]+=this_time[t]
    return [t/n for t in sum_time]


        