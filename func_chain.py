



from sage.all import is_odd,is_prime,sqrt,parent

from func_elliptic import Translate_by_M,Mont_pt_to_lv2,Mont_4torsion_to_lv2null,Pt_to_1m1,Lv2tnp_to_Mont_coeff,Lv2tc_to_Mont_coord,Mont_to_elliptic,Mont_coord_to_point,Decomp_degree
from class_theta_dim1 import Dim1_theta,Dim1_theta_null
from func_proposed_isogeny import Codomain_dim2,Evaluation_dim2_special
from func_elliptic import Is_Elliptic_product
from func_isogeny import CodOne,EvalOne,EvalSq,Product_power_lambda
from class_theta_dim2 import NullCoord,To_null,Coord

import time




def Theta_on_elliptic(M,x,lv2tnp):
    """ 
    M1,M2 are the transformation matrices from Montgomery to Montgomery.
    P1,P2 are any points on Montgomery curve E1,E2.
    tnp_E1,tnp_E2 are theta-null points of E1,E2.
    return the theta coordinate of P=(P1,P2) on E1*E2.
    """
    mon_x = Translate_by_M(M,x)
    lv2_x = Mont_pt_to_lv2(lv2tnp,mon_x)
    return Dim1_theta(lv2_x,1,parent(x[0]))





def Preparation_dim1_theta(E1,E2,P1,P2,Q1,Q2,x1_list:list,N:int):
    assert (P1 in E1)
    assert (P2 in E2)
    assert (Q1 in E1)
    assert (Q2 in E2)
    assert (P1.order()==N)
    assert (P2.order()==N)
    assert (Q1.order()==N)
    assert (Q2.order()==N)
    assert (is_odd(N))
    l=Decomp_degree(N)[0]
    k=N//l
    assert (is_prime(l))
    assert (is_odd(l))
    # on E_1
    e1=k*P1
    f1=k*Q1
    S1d,T1d=E1.torsion_basis(4)
    M1=Pt_to_1m1(T1d)    # M_1: Q_1d|->(1:-1)
    mon_P1d  = Translate_by_M(M1,S1d)
    tnp_E1  = Mont_4torsion_to_lv2null(mon_P1d)
    lv2_e1   = Theta_on_elliptic(M1,e1   ,tnp_E1)
    lv2_f1   = Theta_on_elliptic(M1,f1   ,tnp_E1)
    lv2_e1f1 = Theta_on_elliptic(M1,e1+f1,tnp_E1)
    ext_P1   =[Theta_on_elliptic(M1,P1   ,tnp_E1),
               Theta_on_elliptic(M1,P1+e1,tnp_E1),
               Theta_on_elliptic(M1,P1+f1,tnp_E1)]
    ext_Q1   =[Theta_on_elliptic(M1,Q1   ,tnp_E1),
               Theta_on_elliptic(M1,Q1+e1,tnp_E1),
               Theta_on_elliptic(M1,Q1+f1,tnp_E1)]
    first_basis_E1 =[lv2_e1,lv2_f1,lv2_e1f1]
    ext_x_E1_list =[[Theta_on_elliptic(M1,x1   ,tnp_E1),
                      Theta_on_elliptic(M1,x1+e1,tnp_E1),
                      Theta_on_elliptic(M1,x1+f1,tnp_E1)] for x1 in x1_list]
    # on E_2    
    e2=k*P2
    f2=k*Q2
    S2d,T2d=E2.torsion_basis(4)
    M2=Pt_to_1m1(T2d)    # M_1: Q_1d|->(1:-1)
    mon_P2d      = Translate_by_M(M2,S2d)
    tnp_E2      = Mont_4torsion_to_lv2null(mon_P2d)
    lv2_e2   = Theta_on_elliptic(M2,e2   ,tnp_E2)
    lv2_f2   = Theta_on_elliptic(M2,f2   ,tnp_E2)
    lv2_e2f2 = Theta_on_elliptic(M2,e2+f2,tnp_E2)
    ext_P2   =[Theta_on_elliptic(M2,P2   ,tnp_E2),
               Theta_on_elliptic(M2,P2+e2,tnp_E2),
               Theta_on_elliptic(M2,P2+f2,tnp_E2)]
    ext_Q2   =[Theta_on_elliptic(M2,Q2   ,tnp_E2),
               Theta_on_elliptic(M2,Q2+e2,tnp_E2),
               Theta_on_elliptic(M2,Q2+f2,tnp_E2)]
    first_basis_E2 =[lv2_e2,lv2_f2,lv2_e2f2]
    return tnp_E1,first_basis_E1,ext_P1,ext_Q1,ext_x_E1_list,tnp_E2,first_basis_E2,ext_P2,ext_Q2
           







def Odd_deg_isogeny_chain(
    tnp_E1:Dim1_theta_null,first_basis_E1:list,ext_P1:list,ext_Q1:list,ext_x_E1_list:list,
    tnp_E2:Dim1_theta_null,first_basis_E2:list,ext_P2:list,ext_Q2:list,
    N:int,One_or_Square="One",proposed_or_exsisited="Proposed"):
    """ 
    function for attack (2/3): Main!!
    Compute (N_A,N_A)-isogeny by splitting to isogeny chain.
    At the last step, we don't split to the product of elliptic curves in this function.
    """
    assert (One_or_Square   in {"One"    ,"Square"})
    assert (proposed_or_exsisited in {"Proposed","Existing"})
    fac=Decomp_degree(N)
    print("isogeny chain:",fac)
    #################################################
    #1st step.
    l=fac[0]
    print("ell=",l)
    each_strat=time.perf_counter()
    #Codomain.---------------------------------------   
    tc_0,_,exc_h_lincom_lsq_E2,cod_coeff_E1,cod_coeff_E2=Codomain_dim2(
        tnp_E1,tnp_E2,first_basis_E1,first_basis_E2,l,"One")    
    assert not Is_Elliptic_product(tc_0)[0]
    assert(len(dict(cod_coeff_E1))==l**2)
    assert(len(dict(cod_coeff_E2))==l**2)
    # Evaluation of kernel. ---------------------------------------
    
    # tc_f1=Evaluation_dim2_general(
    #             tnp_E1,first_basis_E1,ext_P1,cod_coeff_E1,
    #             tnp_E2,first_basis_E2,ext_P2,cod_coeff_E2,l,One_or_Square)    
    # tc_f2=Evaluation_dim2_general(
    #             tnp_E1,first_basis_E1,ext_Q1,cod_coeff_E1,
    #             tnp_E2,first_basis_E2,ext_Q2,cod_coeff_E2,l,One_or_Square)
    
    tnp_E1E2=To_null(tnp_E1.Product_theta(tnp_E2))
    first_basis_E1E2=[
        first_basis_E1[i].Product_theta(first_basis_E2[i]) for i in range(0,3)]
    for i in range(0,3):
        first_basis_E1E2[i].order=l
    first_basis_E1E2[0].lmd_lpow_value=[
        cod_coeff_E1[(1,0)][u]*cod_coeff_E2[(1,0)][u] for u in range(0,2)]
    first_basis_E1E2[1].lmd_lpow_value=[
        cod_coeff_E1[(0,1)][u]*cod_coeff_E2[(0,1)][u] for u in range(0,2)]
    first_basis_E1E2[2].lmd_lpow_value=[
        cod_coeff_E1[(1,1)][u]*cod_coeff_E2[(1,1)][u] for u in range(0,2)]
    ext_P1P2=[ext_P1[i].Product_theta(ext_P2[i]) for i in range(0,3)]
    ext_Q1Q2=[ext_Q1[i].Product_theta(ext_Q2[i]) for i in range(0,3)]
    lmd_data=Product_power_lambda(first_basis_E1E2)
    tc_f1=EvalOne(tnp_E1E2,first_basis_E1E2,ext_P1P2,lmd_data)
    tc_f2=EvalOne(tnp_E1E2,first_basis_E1E2,ext_Q1Q2,lmd_data)   
    # assert tc_0.Is_order(tc_f1,N//l)
    # assert tc_0.Is_order(tc_f2,N//l)
    # assert tc_0.Is_on_Kummer_eq(tc_f1)
    # assert tc_0.Is_on_Kummer_eq(tc_f2)
    
    # Main part !!!!! ----------------------------------
    # Evaluation of (x,0).
    if proposed_or_exsisited=="Proposed":
        tc_fx_list =[Evaluation_dim2_special(
                    tnp_E1,first_basis_E1,ext_x_E1,cod_coeff_E1,exc_h_lincom_lsq_E2,l,One_or_Square)
                    for ext_x_E1 in ext_x_E1_list]
    else: #geenral(=exisited method).
        ext_0_E2=[tnp_E2,first_basis_E2[0],first_basis_E2[1]]
        ext_x0_list=[[
            ext_x_E1[i].Product_theta(ext_0_E2[i]) for i in range(0,3)]
                     for ext_x_E1 in ext_x_E1_list]
        if One_or_Square=="One":
            tc_fx_list=[EvalOne(tnp_E1E2,first_basis_E1E2,ext_x0,lmd_data)
                    for ext_x0 in ext_x0_list]
        else: #Square
            tc_fx_list=[EvalSq(tnp_E1E2,first_basis_E1E2,ext_x0,lmd_data)
                    for ext_x0 in ext_x0_list]
    # for tc_fx in tc_fx_list:
    #     assert tc_0.Is_on_Kummer_eq(tc_fx)
    each_end=time.perf_counter()
    print("Time:",each_end-each_strat,"sec")
    #################################################
    #2nd~last.
    s=fac[0]
    for i in range(1,len(fac)):
        each_strat=time.perf_counter()
        l=fac[i]
        assert(is_prime(l))
        k=N//(s*l)
        assert(s*l*k==N)
        print("ell=",l)
        # assert tc_0.Is_order(tc_f1,l*k)
        # assert tc_0.Is_order(tc_f2,l*k)
        # assert not tc_f1.Peq(tc_f2)
        tc_f12=tc_0.Normal_Add(tc_f1,tc_f2,1)
        assert(type(tc_f12)==Coord)
        # assert tc_0.Is_order(tc_f12,l*k)
        #construct kernel e_1,e_2,e_1+e_2.
        tc_e1 =tc_0.Mult(tc_f1 ,k)
        tc_e2 =tc_0.Mult(tc_f2 ,k)
        tc_e12=tc_0.Mult(tc_f12,k)
        tc_e1.order=l
        tc_e2.order=l
        tc_e12.order=l
        # assert tc_0.Is_order(tc_e1,l)
        # assert tc_0.Is_order(tc_e2,l)
        # assert tc_0.Is_order(tc_e12,l)
        tc_f1pe1=tc_0.Mult(tc_f1,k+1)               #f_1+e_1
        tc_f1pe2=tc_0.Kxpy_xpy(k,tc_f2,tc_f1,tc_f12)#f_1+e_2
        tc_f2pe1=tc_0.Kxpy_xpy(k,tc_f1,tc_f2,tc_f12)#f_2+e_1
        tc_f2pe2=tc_0.Mult(tc_f2,k+1)               #f_2+e_2
        ext_x_list=[]
        for i in range(0,len(tc_fx_list)):
            tc_x=tc_fx_list[i]
            tc_xpe1 =tc_0.Normal_Add(tc_x,tc_e1,1)#x+e_1
            tc_xpe2 =tc_0.Compatible_Add(tc_x,tc_e2,tc_xpe1,tc_e12)#x+e_2
            ext_x_list.append([tc_x,tc_xpe1,tc_xpe2])
        basis  =[tc_e1,tc_e2   ,tc_e12  ]
        f1_list=[tc_f1,tc_f1pe1,tc_f1pe2]
        f2_list=[tc_f2,tc_f2pe1,tc_f2pe2]
        #---------------------------------------------------------------
        tc_cd0=CodOne(tc_0,basis)
        lmd_data=Product_power_lambda(basis)
        tc_f1=EvalOne(tc_0,basis,f1_list,lmd_data)
        tc_f2=EvalOne(tc_0,basis,f2_list,lmd_data)
        tc_fx_list=[EvalOne(tc_0,basis,ext_x_list[i],lmd_data) for i in range(0,len(tc_fx_list))]
        tc_0 =tc_cd0
        s*=l
        each_end=time.perf_counter()
        #print("Time:",each_end-each_strat,"sec")
    #=====================================
    assert(Is_Elliptic_product(tc_0)[0])
    assert(tc_0.Is_same_proj(tc_f1))
    assert(tc_0.Is_same_proj(tc_f1))
    return tc_0,tc_fx_list







def Theta_Split(lv2tnp:NullCoord,lv2_x_list:list,zeta_4):
    """ 
    for given theta null point on product of elliptic curves, 
    then this function gives theta null of product of theta.
    cf[DMPR23]p18.
    """
    assert(zeta_4**2==-1)
    assert(Is_Elliptic_product(lv2tnp)[0])
    if Is_Elliptic_product(lv2tnp)[1]:  #already split.
        return lv2tnp,lv2_x_list
    else:
        _,_,i,j=Is_Elliptic_product(lv2tnp)
    if (i,j)==(0,0):
        lv2tnp.numer[2]*=zeta_4
        lv2tnp.numer[3]*=zeta_4
        for k in range(0,len(lv2_x_list)):
            lv2_x_list[k].numer[2] *=zeta_4
            lv2_x_list[k].numer[3] *=zeta_4
    if Is_Elliptic_product(lv2tnp)[1]:
        return lv2tnp,lv2_x_list
    else:
        _,_,i,j=Is_Elliptic_product(lv2tnp) #(i,j) is zero even theta.
    if i!=0 and j==0:
        assert(i!=0)
        lv2tnp=To_null(lv2tnp.Hadamard())
        for k in range(0,len(lv2_x_list)):
            lv2_x_list[k] = lv2_x_list[k].Hadamard()
    if Is_Elliptic_product(lv2tnp)[1]:
        return lv2tnp,lv2_x_list
    else:
        _,_,i,j=Is_Elliptic_product(lv2tnp)#(i,j) is zero even theta.
    assert(j!=0)
    if j==1:
        (lv2tnp.numer[1],lv2tnp.numer[3])=(lv2tnp.numer[3],lv2tnp.numer[1])
        for k in range(0,len(lv2_x_list)):
            (lv2_x_list[k].numer[1] ,lv2_x_list[k].numer[3])=(
                lv2_x_list[k].numer[3] ,lv2_x_list[k].numer[1])
    elif j==2:
        (lv2tnp.numer[2],lv2tnp.numer[3])=(lv2tnp.numer[3],lv2tnp.numer[2])
        for k in range(0,len(lv2_x_list)):
            (lv2_x_list[k].numer[2] ,lv2_x_list[k].numer[3])=(
                lv2_x_list[k].numer[3] ,lv2_x_list[k].numer[2])
    if Is_Elliptic_product(lv2tnp)[1]:
        return lv2tnp,lv2_x_list
    else:
        _,_,i,j=Is_Elliptic_product(lv2tnp)#(i,j) is zero even theta.
    assert(i==0 and j==3)
    lv2tnp.numer[1]*=zeta_4
    lv2tnp.numer[2]*=zeta_4
    for k in range(0,len(lv2_x_list)):
        lv2_x_list[k].numer[1]*=zeta_4
        lv2_x_list[k].numer[2]*=zeta_4
    assert(Is_Elliptic_product(lv2tnp)[1])
    return lv2tnp,lv2_x_list
    



def Splitting(lv2tnp:NullCoord,tc_fx_list:list):
    """ 
    function for attack (3/3)
    From  the last step of isongey chain, 
    we get the secret key of Bob in B-SIDH.
    """
    base_field=lv2tnp.field
    zeta_4=sqrt(base_field(-1))
    lv2tnp,tc_fx_list=Theta_Split(lv2tnp,tc_fx_list,zeta_4) 
    assert(Is_Elliptic_product(lv2tnp)[0])
    assert(Is_Elliptic_product(lv2tnp)[1])
    assert(lv2tnp.numer[0]*lv2tnp.numer[3]==lv2tnp.numer[1]*lv2tnp.numer[2])
    #split to elliptic curves.----------------------
    lv2tnp_Ecd_1=Dim1_theta_null([lv2tnp.numer[0],lv2tnp.numer[1]],1,base_field)
    lv2tnp_Ecd_2=Dim1_theta_null([lv2tnp.numer[0],lv2tnp.numer[2]],1,base_field)
    tc_fx_E1_list=[]
    tc_fx_E2_list=[]
    for tc_fx in tc_fx_list:
        tc_fx_E1_list.append(Dim1_theta([
            tc_fx.numer[0],tc_fx.numer[1]],1,base_field))
        tc_fx_E2_list.append(Dim1_theta([
            tc_fx.numer[0],tc_fx.numer[2]],1,base_field))
    #change to Montgomery.----------------------
    A1_cod=Lv2tnp_to_Mont_coeff(lv2tnp_Ecd_1)
    A2_cod=Lv2tnp_to_Mont_coeff(lv2tnp_Ecd_2)
    mont_fx_E1_list=[]
    mont_fx_E2_list=[]
    for k in range(0,len(tc_fx_list)):
        mont_fx_E1_list.append(Lv2tc_to_Mont_coord(lv2tnp_Ecd_1,tc_fx_E1_list[k]))
        mont_fx_E2_list.append(Lv2tc_to_Mont_coord(lv2tnp_Ecd_2,tc_fx_E2_list[k]))
    #Montgomery to Sage's elliptic.----------------------
    E1_cod=Mont_to_elliptic(A1_cod)
    E2_cod=Mont_to_elliptic(A2_cod)
    pt_fx_E1_list=[]
    pt_fx_E2_list=[]
    for k in range(0,len(tc_fx_list)):
        pt_fx_E1_list.append(Mont_coord_to_point(E1_cod,A1_cod,mont_fx_E1_list[k]))
        pt_fx_E2_list.append(Mont_coord_to_point(E2_cod,A2_cod,mont_fx_E2_list[k]))
    return E1_cod,E2_cod,pt_fx_E1_list,pt_fx_E2_list
    




def Total_computation(E1,E2,P1,P2,Q1,Q2,fx_list:list,D:int,One_or_Square,proposed_or_exsisited):
    """ 
    This function is the main function that integrates all the individual steps.
    """
    time_strat=time.perf_counter()
    #compute theta coordinate to start isogeny chain.
    tnp_E1,first_basis_E1,ext_P1,ext_Q1,ext_x_E1_list,tnp_E2,first_basis_E2,ext_P2,ext_Q2=Preparation_dim1_theta(E1,E2,P1,P2,Q1,Q2,fx_list,D)
    #Main: isogeny chain.
    tnp_cod,tc_x_cod_list=Odd_deg_isogeny_chain(tnp_E1,first_basis_E1,ext_P1,ext_Q1,ext_x_E1_list,tnp_E2,first_basis_E2,ext_P2,ext_Q2,D,One_or_Square,proposed_or_exsisited)
    #split to elliptic curves.
    result=Splitting(tnp_cod,tc_x_cod_list)
    time_end=time.perf_counter()
    print("Total time:",time_end-time_strat,"sec")
    return result


