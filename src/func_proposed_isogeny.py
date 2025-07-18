
#################################################################
# In this module, we give some functions for (ell,ell)-isogeny calculation where ell is odd prime number.
# There are 3 types functions: Codomain, Evaluation(general), Evaluation(special). 
###################################################################

from sage.all import prod

from func_lin_combi   import Set_H_ell,Half_LinCom_dim1,XpLinCom_dim1,Comp_H_ell,Remain_Half_coeff_without0
from func_fraction    import Dict_common_denom_len3,Multpower_straight,Multpower_sq,Common_denom_frac
from class_theta_dim1 import Dim1_theta,Dim1_theta_null
from class_theta_dim2 import Coord,NullCoord
from func_existing_isogeny import Sum_of_square


        



def Prod_nn_nn(nnd1:list,nnd2:list):
    """
    "nnd" is 1-dimensional theta coordinate holded as [numerator at 0,numerator at 1,(common) denominator].
    For nnd1, nnd2, we give the product theta coordinate in dim 2.
    """
    assert len(nnd1)==3
    assert len(nnd2)==3
    prod=[nnd1[0]*nnd2[0],
          nnd1[1]*nnd2[0],
          nnd1[0]*nnd2[1],
          nnd1[1]*nnd2[1]]
    return prod
    
#############################################################################
#############################################################################


def Lmd_ell_sq_dim1(tc_ld_e:Dim1_theta,tc_ldp1_e:Dim1_theta):
    """ compute lmd^{\l} to use normalization of affine lift of basis."""
    assert (tc_ld_e.Peq(tc_ldp1_e))
    lmd_lsq_1=[tc_ld_e.th0*tc_ldp1_e.denom,tc_ld_e.denom*tc_ldp1_e.th0]
    #lmd_lsq_2=[tc_ld_e.th1*tc_ldp1_e.denom,tc_ld_e.denom*tc_ldp1_e.th1]
    #assert lmd_lsq_1[0]/lmd_lsq_1[1]==lmd_lsq_2[0]/lmd_lsq_2[1]
    return lmd_lsq_1




def Coeff_term_cod_dim1(
        ld_e1 :Dim1_theta,ldp1_e1 :Dim1_theta,
        ld_e2 :Dim1_theta,ldp1_e2 :Dim1_theta,
        ld_e12:Dim1_theta,ldp1_e12:Dim1_theta,l:int):
    """ compute coefficients to use the normalization for codomain calculation in dimension 1."""
    #pick up.
    lmd1_lpow =Lmd_ell_sq_dim1(ld_e1 ,ldp1_e1)
    lmd2_lpow =Lmd_ell_sq_dim1(ld_e2 ,ldp1_e2)
    lmd12_lpow=Lmd_ell_sq_dim1(ld_e12,ldp1_e12)
    lmd_div_lpow=[
        lmd12_lpow[0]*lmd1_lpow[1]*lmd2_lpow[1],
        lmd12_lpow[1]*lmd1_lpow[0]*lmd2_lpow[0]]
    [lmd1_lpow,lmd2_lpow,lmd_div_lpow],den=Common_denom_frac(
        [lmd1_lpow,lmd2_lpow,lmd_div_lpow])
    return lmd1_lpow,lmd2_lpow,lmd_div_lpow,den





def Excellent_term_CodOne_dim1(tc_0:Dim1_theta_null,ext_basis:list,l:int):
    """ normalization in codomain in dimension 1."""
    #linear combination.
    h_lincom=Half_LinCom_dim1(tc_0,ext_basis,l)
    h_lincom_lsq=dict()
    for (k1,k2) in Set_H_ell(l)|{(0,0)}:
        h_lincom_lsq[(k1,k2)]=[h_lincom[(k1,k2)].th0**l,
                               h_lincom[(k1,k2)].th1**l, 
                               h_lincom[(k1,k2)].denom**l]
    #coefficients.
    ld=(l-1)//2
    assert (h_lincom[(ld,0)].Peq(h_lincom[(ld+1,0)]))
    lmd1_lpow,lmd2_lpow,lmd_div_lpow,den=Coeff_term_cod_dim1(
                                h_lincom[(ld,0)] ,h_lincom[(ld+1,0)],
                                h_lincom[(0,ld)] ,h_lincom[(0,ld+1)],
                                h_lincom[(ld,ld)],h_lincom[(ld+1,ld+1)],l)
    assert(lmd1_lpow[1]==lmd2_lpow[1]==lmd_div_lpow[1]==den)
    den_pow         =Multpower_straight(den,3*(l-1)**2)
    lmd1_lpow_pow   =Multpower_sq(lmd1_lpow[0],l)
    lmd2_lpow_pow   =Multpower_sq(lmd2_lpow[0],l)
    lmd_div_lpow_pow=Multpower_straight(lmd_div_lpow[0],(l-1)**2)
    coeff=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            coeff[(k1,k2)]=[
                lmd1_lpow_pow[k1**2]*lmd2_lpow_pow[k2**2]*lmd_div_lpow_pow[k1*k2],
                den_pow[k1**2+k2**2+k1*k2]]
    #normalization.
    excellent_h_lincom_lsq=dict()
    for (k1,k2) in Set_H_ell(l)|{(0,0)}:
        excellent_h_lincom_lsq[(k1,k2)]=[
            h_lincom_lsq[(k1,k2)][0]*coeff[(k1,k2)][0],
            h_lincom_lsq[(k1,k2)][1]*coeff[(k1,k2)][0],
            h_lincom_lsq[(k1,k2)][2]*coeff[(k1,k2)][1]] 
    excellent_h_lincom_lsq,_=Dict_common_denom_len3(excellent_h_lincom_lsq)
    return excellent_h_lincom_lsq,coeff






def Excellent_term_CodSq_dim1(tc_0:Dim1_theta_null,ext_basis:list,l:int):
    """ normalization in codomain in dimension 1."""
    #take l-th power.
    h_lincom=Half_LinCom_dim1(tc_0,ext_basis,l)
    ld=(l-1)//2
    lmd1_lpow,lmd2_lpow,lmd_div_lpow,den=Coeff_term_cod_dim1(
                                h_lincom[(ld,0)] ,h_lincom[(ld+1,0)],
                                h_lincom[(0,ld)] ,h_lincom[(0,ld+1)],
                                h_lincom[(ld,ld)],h_lincom[(ld+1,ld+1)],l)
    a_u=Sum_of_square(l)
    r=len(a_u)
    ss1=[(a_u[u]*(l-1))%l         for u in range(0,r)]
    tt1=[(a_u[u]*(l-1)-ss1[u])//l for u in range(0,r)]
    max_exp=(l-1)**2+l*(sum([tt1[u]**2 for u in range(0,r)]))-2*(l-1)*(sum([a_u[u]*tt1[u] for u in range(0,r)]))
    assert(max_exp<=r*(l-1))
    assert(lmd1_lpow[1]==lmd2_lpow[1]==lmd_div_lpow[1]==den)
    lmd1_lpow_pow   =Multpower_straight(lmd1_lpow[0]   ,max_exp)
    lmd2_lpow_pow   =Multpower_straight(lmd2_lpow[0]   ,max_exp)
    lmd_div_lpow_pow=Multpower_straight(lmd_div_lpow[0],max_exp)
    den_pow         =Multpower_straight(den            ,max_exp*3)
    #extension of linear combination to all 0<=k1,k2<l.
    lincom=dict()
    lincom[(0,0)]=tc_0
    for (k1,k2) in Set_H_ell(l):
        lincom[(k1,k2)]=h_lincom[(k1,k2)]
    for (k1,k2) in Remain_Half_coeff_without0(l):
        assert(not ((k1,k2) in lincom))
        a1  =l-2*k1
        a2  =l-2*k2
        adiv=l-k1-k2
        aden=a1+a2+adiv
        assert(aden==3*adiv)
        if k1==0:
            assert(not (0,k2) in lincom)
            assert(a2<=0)
            lincom[(0,k2)]=lincom[(0,l-k2)].Mult_frac(
                [den_pow[-a2],lmd2_lpow_pow[-a2]])
        elif k2==0:
            assert(not (k1,0) in lincom)
            assert(a1<=0)
            lincom[(k1,0)]=lincom[(l-k1,0)].Mult_frac(
                [den_pow[-a1],lmd1_lpow_pow[-a1]])
        elif (a1>=0) and (a2>=0):
            assert(not (k1,k2) in lincom)
            lincom[(k1,k2)]=lincom[(l-k1,l-k2)].Mult_frac(
                [lmd1_lpow_pow[a1]*lmd2_lpow_pow[a2]*den_pow[-aden],lmd_div_lpow_pow[-adiv]])
        elif (a1>=0) and (a2<0):
            assert(not (k1,k2) in lincom)
            lincom[(k1,k2)]=lincom[(l-k1,l-k2)].Mult_frac(
                [lmd1_lpow_pow[a1]*den_pow[-aden],lmd2_lpow_pow[-a2]*lmd_div_lpow_pow[-adiv]])
        elif (a1<0) and (a2>=0):
            assert(not (k1,k2) in lincom)
            lincom[(k1,k2)]=lincom[(l-k1,l-k2)].Mult_frac(
                [lmd2_lpow_pow[a2]*den_pow[-aden],lmd1_lpow_pow[-a1]*lmd_div_lpow_pow[-adiv]])
        else:
            assert((a1<0) and (a2<0))
            assert(not (k1,k2) in lincom)
            lincom[(k1,k2)]=lincom[(l-k1,l-k2)].Mult_frac(
                [den_pow[-aden],lmd1_lpow_pow[-a1]*lmd2_lpow_pow[-a2]*lmd_div_lpow_pow[-adiv]])
    #assert(len(lincom)==l**2)
    #compute excellent term.
    excellent_h_lincom_lsq=dict()
    coeff=dict()
    excellent_h_lincom_lsq[(0,0)]=[
        prod([tc_0.numer[(a_u[u]%2)*0] for u in range(0,r)]),
        prod([tc_0.numer[(a_u[u]%2)*1] for u in range(0,r)]),
        (tc_0.denom)**r
    ]
    for (k1,k2) in Set_H_ell(l)|{(0,0)}:
        s1=[(a_u[u]*k1)%l        for u in range(0,r)]
        t1=[(a_u[u]*k1-s1[u])//l for u in range(0,r)]
        s2=[(a_u[u]*k2)%l        for u in range(0,r)]
        t2=[(a_u[u]*k2-s2[u])//l for u in range(0,r)]
        for u in range(0,r):
            assert(a_u[u]*k1==t1[u]*l+s1[u])
            assert(a_u[u]*k2==t2[u]*l+s2[u])
            assert(0<=s1[u]<l)
            assert(0<=s2[u]<l)
        h_1  =k1**2+l*(sum([t1[u]**2    for u in range(0,r)]))-2*k1*(sum([a_u[u]*t1[u] for u in range(0,r)]))
        h_2  =k2**2+l*(sum([t2[u]**2    for u in range(0,r)]))-2*k2*(sum([a_u[u]*t2[u] for u in range(0,r)]))
        h_div=k1*k2+l*(sum([t1[u]*t2[u] for u in range(0,r)]))-k1*sum([a_u[u]*t2[u] for u in range(0,r)])-k2*sum([a_u[u]*t1[u] for u in range(0,r)])
        h_den=h_1+h_2+h_div
        coeff[(k1,k2)]=[lmd1_lpow_pow[h_1]*lmd2_lpow_pow[h_2]*lmd_div_lpow_pow[h_div],den_pow[h_den]]
        excellent_h_lincom_lsq[(k1,k2)]=[
            coeff[(k1,k2)][0]*prod(
                [lincom[(s1[u],s2[u])].numer[(a_u[u]%2)*0] for u in range(0,r)]),
            coeff[(k1,k2)][0]*prod(
                [lincom[(s1[u],s2[u])].numer[(a_u[u]%2)*1] for u in range(0,r)]),
            coeff[(k1,k2)][1]*prod(
                [lincom[(s1[u],s2[u])].denom               for u in range(0,r)])
        ]
    excellent_h_lincom_lsq,_=Dict_common_denom_len3(excellent_h_lincom_lsq)
    #prepare coefficients for evaluation.
    den_pow         =Multpower_straight(den,3*(l-1)**2)
    lmd1_lpow_pow   =Multpower_sq(lmd1_lpow[0],l)
    lmd2_lpow_pow   =Multpower_sq(lmd2_lpow[0],l)
    lmd_div_lpow_pow=Multpower_straight(lmd_div_lpow[0],(l-1)**2)
    coeff_for_eval=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            coeff_for_eval[(k1,k2)]=[
                lmd1_lpow_pow[k1**2]*lmd2_lpow_pow[k2**2]*lmd_div_lpow_pow[k1*k2],
                den_pow[k1**2+k2**2+k1*k2]]
    return excellent_h_lincom_lsq,coeff_for_eval





def Codomain_dim2(tc_0_E1:Dim1_theta_null,tc_0_E2:Dim1_theta_null,ext_basis_E1:list,ext_basis_E2:list,l:int,method="One"):
    """ 
    Compute theta-null point of the codomain for (l,l)-isongey.
    Here, e1,e2 is a basis of the kernel in (E1*E2)[l].
    tc_0_E1 is the theta-null point of E1 and ext_basis_E1 is the first component of e1,e2,e1+e2.
    """
    assert(method in {"One","Square"})
    # X,_=Excellent_term_CodOne_dim1(tc_0_E1,ext_basis_E1,l)
    # Y,_=Excellent_term_CodSq_dim1(tc_0_E1,ext_basis_E1,l)
    # print(X[(0,0)][0]/X[(0,0)][1]==Y[(0,0)][0]/Y[(0,0)][1])
    if method=="One":
        exc_h_lincom_lsq_E1,cod_coeff_E1=Excellent_term_CodOne_dim1(
            tc_0_E1,ext_basis_E1,l)
        exc_h_lincom_lsq_E2,cod_coeff_E2=Excellent_term_CodOne_dim1(
            tc_0_E2,ext_basis_E2,l)
    else:
        exc_h_lincom_lsq_E1,cod_coeff_E1=Excellent_term_CodSq_dim1(
            tc_0_E1,ext_basis_E1,l)
        exc_h_lincom_lsq_E2,cod_coeff_E2=Excellent_term_CodSq_dim1(
            tc_0_E2,ext_basis_E2,l)
    #assert(len(dict(cod_coeff_E1))==l**2)
    #assert(len(dict(cod_coeff_E2))==l**2)
    tc_f0=Prod_nn_nn(exc_h_lincom_lsq_E1[(0,0)],exc_h_lincom_lsq_E2[(0,0)])
    for (k1,k2) in Set_H_ell(l):
        sum_part=Prod_nn_nn(
            exc_h_lincom_lsq_E1[(k1,k2)],
            exc_h_lincom_lsq_E2[(k1,k2)])
        for i in range(0,4):
            tc_f0[i]+=2*sum_part[i]
    return NullCoord(tc_f0,1,tc_0_E1.field),exc_h_lincom_lsq_E1,exc_h_lincom_lsq_E2,cod_coeff_E1,cod_coeff_E2




######################################################################
######################################################################
######################################################################
######################################################################




def Product_power_lambda_dim1(lmd1_lpow,lmd2_lpow,lmd12_lpow,l:int):
    lmd_div_lpow=[lmd12_lpow[0]*lmd1_lpow[1]*lmd2_lpow[1],
                  lmd12_lpow[1]*lmd1_lpow[0]*lmd2_lpow[0]]
    [lmd1_lpow,lmd2_lpow,lmd_div_lpow],den=Common_denom_frac(
        [lmd1_lpow,
         lmd2_lpow,lmd_div_lpow])
    den_pow         =Multpower_straight(den            ,3*(l-1)**2)
    lmd1_lpow_pow   =Multpower_straight(lmd1_lpow[0]   ,l**2)
    lmd2_lpow_pow   =Multpower_straight(lmd2_lpow[0]   ,l**2)
    lmd_div_lpow_pow=Multpower_straight(lmd_div_lpow[0],(l-1)**2)
    lmd_pow_product=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            lmd_pow_product[(k1,k2)]=[
                lmd1_lpow_pow[k1**2]*lmd2_lpow_pow[k2**2]*lmd_div_lpow_pow[k1*k2],
                den_pow[k1**2+k2**2+k1*k2]]
    lm_1_lsq=[lmd1_lpow_pow[l],den_pow[l]]
    lm_2_lsq=[lmd2_lpow_pow[l],den_pow[l]]
    #assert(len(lmd_pow_product)==l**2)
    return [lmd_pow_product,lm_1_lsq,lm_2_lsq]




def MudivLmd_ell_sq_dim1(tc_x:Dim1_theta,tc_xple:Dim1_theta,lmd_lsq,l):
    """ compute (mu/lmd)^l."""
    #(mu/lmd)^l*(lmd^l)^l.
    mdll=[tc_x.th0*tc_xple.denom,tc_x.denom*tc_xple.th0]
    #(mu/lmd)^l.
    mdl=[mdll[0]*lmd_lsq[1]**l,mdll[1]*lmd_lsq[0]**l]
    return mdl




def Coeff_term_eval_dim1(cod_coeff:dict,tc_x:Dim1_theta,tc_xple1:Dim1_theta,tc_xple2:Dim1_theta,l:int):
    """ 
    coefficient appearing in the formula to compute evaluation. 
    we reuse "coeff" we computed in codomain.
    """
    lmd1_lsq=cod_coeff[(1,0)]
    lmd2_lsq=cod_coeff[(0,1)]
    mu1_div_lmd1=MudivLmd_ell_sq_dim1(tc_x,tc_xple1,lmd1_lsq,l)
    mu2_div_lmd2=MudivLmd_ell_sq_dim1(tc_x,tc_xple2,lmd2_lsq,l)
    [mu1_div_lmd1,mu2_div_lmd2],_=Common_denom_frac(
        [mu1_div_lmd1,mu2_div_lmd2])
    assert mu1_div_lmd1[1]==mu2_div_lmd2[1]
    return mu1_div_lmd1,mu2_div_lmd2




def Excellent_term_EvalOne_dim1(tc_0:Dim1_theta_null,ext_basis:list,ext_x:list,l:int,cod_coeff:dict):
    """ 
    normalization in evaluation in dimension 2.
    """
    #assert(len(dict(cod_coeff))==l**2)
    xplincom  =XpLinCom_dim1(tc_0,ext_basis,ext_x,l)
    tc_x      =ext_x[0]
    tc_xple1  =xplincom[(l,0)]
    tc_xple2  =xplincom[(0,l)]
    #powers
    mu1_div_lmd1,mu2_div_lmd2=Coeff_term_eval_dim1(
        cod_coeff,tc_x,tc_xple1,tc_xple2,l)
    mu1_div_lmd1_pow=[
        [mu1_div_lmd1[u]**m1 for u in range(0,2)] for m1 in range(0,l)]
    mu2_div_lmd2_pow=[
        [mu2_div_lmd2[u]**m2 for u in range(0,2)] for m2 in range(0,l)]
    eva_coeff=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            eva_coeff[(k1,k2)]=[
                cod_coeff[(k1,k2)][u]*mu1_div_lmd1_pow[k1][u]*mu2_div_lmd2_pow[k2][u] for u in range(0,2)]
    excl_xplincom_lsq=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            coeff=eva_coeff[(k1,k2)]
            excl_xplincom_lsq[(k1,k2)]=[
                coeff[0]*xplincom[(k1,k2)].th[0]**l,
                coeff[0]*xplincom[(k1,k2)].th[1]**l,
                coeff[1]*xplincom[(k1,k2)].th[2]**l]
    excl_xplincom_lsq,_=Dict_common_denom_len3(excl_xplincom_lsq)
    return excl_xplincom_lsq






def Excellent_term_EvalSq_dim1(tc_0: Dim1_theta_null,ext_basis:list,ext_x:list,l:int,cod_coeff:dict):
    """ compute evaluation by l=a_1^2+...+a_r^2."""
    #assert(len(dict(cod_coeff))==l**2)
    [tc_e1,tc_e2,tc_e12]   =ext_basis
    [tc_x,tc_xpe1,tc_xpe2] =ext_x
    #a_u=Sum_of_square_mod24(l)
    a_u=Sum_of_square(l)
    max_au=max(set(a_u))
    r=len(a_u)
    auxplincom_list=[] #a_u*x+linear combination for 0<=u<r.
    mult_e1  =tc_0.Mult_Mult(tc_e1  ,max_au)
    mult_e2  =tc_0.Mult_Mult(tc_e2  ,max_au)
    mult_e12 =tc_0.Mult_Mult(tc_e12 ,max_au)
    mult_x   =tc_0.Mult_Mult(tc_x   ,max_au)
    mult_xpe1=tc_0.Mult_Mult(tc_xpe1,max_au)
    mult_xpe2=tc_0.Mult_Mult(tc_xpe2,max_au)
    for u in range(0,r):
        au=a_u[u]
        au_e1  =mult_e1  [au]
        au_e2  =mult_e2  [au]
        au_e12 =mult_e12 [au]
        au_x   =mult_x   [au]
        au_xpe1=mult_xpe1[au]
        au_xpe2=mult_xpe2[au]
        # au_e1.order =l
        # au_e2.order =l
        # au_e12.order=l
        if u!=0 and au==a_u[u-1]:
            auxplincom=auxplincom_list[u-1]
        else:
            auxplincom=XpLinCom_dim1(
                tc_0,[au_e1,au_e2,au_e12],[au_x,au_xpe1,au_xpe2],l)
        #assert(len(auxplincom.keys())==l**2+2)
        if au==1:
            xplincom=auxplincom
        auxplincom_list.append(auxplincom)
    if all([a_u[u]!=1 for u in range(0,r)]):
        tc_xple1=tc_0.Kxpy_xpy_dim1(l,tc_e1,tc_x,tc_xpe1)
        tc_xple2=tc_0.Kxpy_xpy_dim1(l,tc_e2,tc_x,tc_xpe2)
    else:
        tc_xple1=tc_0.Diff_add_dim1(
            xplincom[(l-1,0)],tc_e1,xplincom[(l-2,0)])
        tc_xple2=tc_0.Diff_add_dim1(
            xplincom[(0,l-1)],tc_e2,xplincom[(0,l-2)])
    mu1_div_lmd1,mu2_div_lmd2=Coeff_term_eval_dim1(
        cod_coeff,tc_x,tc_xple1,tc_xple2,l)
    assert (mu1_div_lmd1[1]==mu2_div_lmd2[1])
    ml1_lpow_pow =[mu1_div_lmd1[0]**m1 for m1 in range(0,l)]
    ml2_lpow_pow =[mu2_div_lmd2[0]**m2 for m2 in range(0,l)]
    mudivlpow_den=[mu2_div_lmd2[1]**m  for m  in range(0,2*l)]
    excl_xplincom_lsq=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            eva_coeff=[
                cod_coeff[(k1,k2)][0]*ml1_lpow_pow[k1]*ml2_lpow_pow[k2]
                ,cod_coeff[(k1,k2)][1]*mudivlpow_den[k1+k2]
                ]
            excl_xplincom_lsq[(k1,k2)]=[
                eva_coeff[0]*prod(
                    [auxplincom_list[u][(k1,k2)].numer[(a_u[u]%2)*0] for u in range(0,r)]),
                eva_coeff[0]*prod(
                    [auxplincom_list[u][(k1,k2)].numer[(a_u[u]%2)*1] for u in range(0,r)]),
                eva_coeff[1]*prod(
                    [auxplincom_list[u][(k1,k2)].denom               for u in range(0,r)])]
    excl_xplincom_lsq,_=Dict_common_denom_len3(excl_xplincom_lsq)
    return excl_xplincom_lsq
    
    




def Evaluation_dim2_general(
    tc_0_E1:Dim1_theta_null,ext_basis_E1:list,ext_x_E1:list,cod_coeff_E1:dict,tc_0_E2:Dim1_theta_null,ext_basis_E2:list,ext_x_E2:list,cod_coeff_E2:dict,l:int,method="One"):
    """ 
    evaluation for the general case, i.e., x=(x_E1,x_E2).
    The input "coeff" is reused value we computed when we calculated codomain.
    """ 
    #assert(len(dict(cod_coeff_E1))==l**2)
    #assert(len(dict(cod_coeff_E2))==l**2)
    if method=="One":
        excl_xplincom_lsq_E1=Excellent_term_EvalOne_dim1(
            tc_0_E1,ext_basis_E1,ext_x_E1,l,cod_coeff_E1)
        excl_xplincom_lsq_E2=Excellent_term_EvalOne_dim1(
            tc_0_E2,ext_basis_E2,ext_x_E2,l,cod_coeff_E2)
    elif method=="Square":
        excl_xplincom_lsq_E1=Excellent_term_EvalSq_dim1(
            tc_0_E1,ext_basis_E1,ext_x_E1,l,cod_coeff_E1)
        excl_xplincom_lsq_E2=Excellent_term_EvalSq_dim1(
            tc_0_E2,ext_basis_E2,ext_x_E2,l,cod_coeff_E2)
    else:
        print("method is not correct.")
        assert(False)
    tc_fx=[0,0,0,0]
    for k1 in range(0,l):
        for k2 in range(0,l):
            sum_part=Prod_nn_nn(
                excl_xplincom_lsq_E1[(k1,k2)],
                excl_xplincom_lsq_E2[(k1,k2)])
            for i in range(0,4):
                tc_fx[i]+=sum_part[i]
    return Coord(tc_fx,1,tc_0_E1.field)





def Evaluation_dim2_special(
    tc_0_E1:Dim1_theta_null,ext_basis_E1:list,ext_x_E1:list,cod_coeff_E1:dict,exc_h_lincom_lsq_E2:dict,l:int,method="One"):
    """ 
    Evaluation for the case of x=(x_E1,0_E2).
    """
    assert (method in {"One","Square"})
    if method=="One":
        excl_xplincom_lsq_E1=Excellent_term_EvalOne_dim1(
            tc_0_E1,ext_basis_E1,ext_x_E1,l,cod_coeff_E1)
    else:
        excl_xplincom_lsq_E1=Excellent_term_EvalSq_dim1(
            tc_0_E1,ext_basis_E1,ext_x_E1,l,cod_coeff_E1)
    excl_xplincom_lsq_E12=dict()
    for (k1,k2) in Set_H_ell(l)|{(0,0)}:
        excl_xplincom_lsq_E12[(k1,k2)]=Prod_nn_nn(
            excl_xplincom_lsq_E1[(k1,k2)],exc_h_lincom_lsq_E2[(k1,k2)])
    c_H_ell=Comp_H_ell(l)
    for (k1,k2) in c_H_ell:
        (k1d,k2d)=c_H_ell[(k1,k2)]
        excl_xplincom_lsq_E12[(k1,k2)]=Prod_nn_nn(
            excl_xplincom_lsq_E1[(k1,k2)],exc_h_lincom_lsq_E2[(k1d,k2d)])
    #assert(len(excl_xplincom_lsq_E12)==l**2)
    tc_fx=[0,0,0,0]
    for k1 in range(0,l):
        for k2 in range(0,l):
            for i in range(0,4):
                tc_fx[i]+=excl_xplincom_lsq_E12[(k1,k2)][i]
    return Coord(tc_fx,1,tc_0_E1.field)

