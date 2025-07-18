

###########################################################
# Functions about elliptic curves.
###########################################################

from sage.all import order,EllipticCurveIsogeny,EllipticCurve,parent,sqrt,factor

from class_theta_dim1 import Dim1_theta_null,Dim1_theta
from class_theta_dim2 import Coord,NullCoord
import itertools


#Cyclic isogeny=========================

def Decomp_degree(N:int):
    """ decompose a given integer to prime factors as list."""
    fac=list(factor(N))
    compo_list=[]
    for i in range(0,len(fac)):
        for counter in range(0,fac[i][1]):
            compo_list.append(fac[i][0])
    return compo_list[::-1] #reverse the order of the list.




def Elliptic_Cyclic(E,ker,P,Q):
    """ 
    Cyclic isogeny from E with kernel "ker".
    Then, compute the image of P,Q.
    """
    deg=order(ker)
    fac_list=Decomp_degree(deg)
    #print("Elliptic isogeny",fac_list)
    s=1
    for l in fac_list:
        k=deg//(s*l)
        decomp_ker=k*ker
        isogeny=EllipticCurveIsogeny(E,decomp_ker)
        E=isogeny.codomain()
        ker=isogeny(ker)
        P=isogeny(P)
        Q=isogeny(Q)
        s*=l
    return E,P,Q






#product theta coordinate =========================================================


def Product_theta(lv2_1:list,lv2_2:list):
    """ compute product theta coordinate."""
    assert(len(lv2_1)==2)
    assert(len(lv2_2)==2)
    lv2=[lv2_1[0]*lv2_2[0],
         lv2_1[0]*lv2_2[1],
         lv2_1[1]*lv2_2[0],
         lv2_1[1]*lv2_2[1]]
    return lv2



def Theta_Hadamard(lv2_x:list):
    assert(len(lv2_x)==4)
    x=lv2_x[0]
    y=lv2_x[1]
    z=lv2_x[2]
    w=lv2_x[3]
    dula_lv2_x=[x+y+z+w,
                x-y+z-w,
                x+y-z-w,
                x-y-z+w]
    return dula_lv2_x





def Product_term(lv2_x,i,chi):
    """ 
    cf.[LR16]section5
    output sum_{t}(chi(t)lv2_x{i+t}lv2_x(t))
    """
    assert(type(lv2_x)==list)
    K=parent(lv2_x[0])
    sum=K(0)
    for t in range(0,4):
        chi_t=(-1)**((chi//2)*(t//2)+(chi%2)*(t%2))
        ipt=2*((i//2+t//2)%2)+(i%2+t%2)%2
        sum+=(chi_t*lv2_x[ipt]*lv2_x[t])
    return sum



def Lv2tnp_to_lv22tnpsq(lv2tnp:NullCoord):
    """ 
    from lv2 theta-null point to lv(2,2) theta null squard.
    """
    assert(type(lv2tnp)==NullCoord)
    lv22tnpsq={}
    for chi in range(0,4):
        for i in range(0,4):
            lv22tnpsq[(chi,i)]=Product_term(lv2tnp.numer,i,chi)
    return lv22tnpsq
            



def Is_Elliptic_product(lv2tnp:NullCoord):
    """ 
    For level 2 theta-null point, check if the base variety (without theta-structure) is a product of elliptic curves.
    """
    assert(type(lv2tnp)==NullCoord)
    even_theta_lv22={(i,j) for i,j in itertools.product(range(0,4),range(0,4)) 
                     if ((i//2)*(j//2)+(i%2)*(j%2))%2==0}
    odd_theta_lv22 ={(i,j) for i,j in itertools.product(range(0,4),range(0,4)) 
                     if ((i//2)*(j//2)+(i%2)*(j%2))%2==1}
    lv22tnpsq=Lv2tnp_to_lv22tnpsq(lv2tnp)
    num_zero=len({k for k in lv22tnpsq if lv22tnpsq[k]==0})
    assert(num_zero==6 or num_zero==7)
    for i_chi in odd_theta_lv22: #odd theta=0
        i=i_chi[0]
        chi=i_chi[1]
        assert(lv22tnpsq[(i,chi)]==0)
    if lv22tnpsq[(3,3)]==0:
        return True,True,0,0 #elliptic product with product theta structure.
    for i_chi in even_theta_lv22: #odd theta=0
        i=i_chi[0]
        chi=i_chi[1]
        if lv22tnpsq[(i,chi)]==0:
            return True,False,i,chi #elliptic product with non-product theta.
    return False,0,0,0 #Jacobian.
   





#Montgomery-to-Montgomery=================================

def Pt_to_1m1(P):
    """ 
    From Montgomery curve with 4-torsion point P,
    output transformation matrix M to other Montgomery curve with 4-torsion point (-1:1).
    """
    assert(2*P!=0)
    assert(4*P==0)
    x=P[0]
    z=P[2]
    u=(2*P)[0]
    w=(2*P)[2]
    M=[-z*w,z*u,0,x*w-u*z]
    return M



def Translate_by_M(M,P:list):
    """ 
    For point P on Montgomery curve, output the image of P under the translation by M.
    """
    x=P[0]
    z=P[2]
    new_x=M[0]*x+M[1]*z
    new_z=M[2]*x+M[3]*z
    return [new_x,new_z]


#Montgomery-to-Thata =================================

def Mont_4torsion_to_lv2null(T_1:list):
    """
    From Montgomery to theta.
    Given 4-torsion point T'_1 s.t. {T'_1;T'_2} is a symeplectic basis with T'_2=(-1:1) in Montgomery coordinate, 
    output theta-null point of the induced theta structure.
    """
    assert(type(T_1)==list)
    assert(len(T_1)==2)
    base_field=parent(T_1[0])
    [r,s]=T_1
    a=r+s
    b=r-s
    lv2tnp=[a,b]
    return Dim1_theta_null(lv2tnp,1,base_field)
    



def Mont_pt_to_lv2(lv2tnp,P:list):
    """ 
    From Montgomery coordinate, output level 2 theta coordinate.
    """
    assert(type(P)==list)
    if P==[0,0]:
        return lv2tnp
    [a,b]=lv2tnp.numer
    [x,z]=P
    lv2_x=[a*(x-z),b*(x+z)]
    return lv2_x



#Theta-to-Montgomery=================================



def Lv2tnp_to_Mont_coeff(lv2tnp:Dim1_theta_null):
    """ 
    From  level 2 theta-null point, output Montgomery coefficient A.
    i.e., E: y^2=x^3+A*x^2+x.
    """
    [a,b]=lv2tnp.numer
    m_coeff=-2*(a**4+b**4)/(a**4-b**4)
    return m_coeff




def Lv2tc_to_Mont_coord(lv2tnp:Dim1_theta_null,lv2_x:Dim1_theta):
    """ 
    From  level 2 theta coordinate, output Montgomery coordinate (x:z).
    """
    [a,b]    =lv2tnp.numer
    [th0,th1]=lv2_x.numer
    return [a*th1+b*th0,a*th1-b*th0]



    
#Montgomery to Sage's elliptic curve=================================


def Mont_to_elliptic(coeff):
    """ 
    For Montgomery curve with coefficient A, 
    i.e., E: y^2=x^3+A*x^2+x,
    output Sage's elliptic curve.
    """
    E=EllipticCurve([0,coeff,0,1,0])
    return E 


def Mont_coord_to_point(E,coeff,xz_coord):
    """ 
    For a point on the Montgomery curve, 
    output the Sage's point on elliptic curve.
    """ 
    x_coord=xz_coord[0]/xz_coord[1]
    y_coord=sqrt(x_coord**3+coeff*x_coord**2+x_coord)
    return E([x_coord,y_coord])


