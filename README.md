# $(\ell,\ell)$-Isogeny from a Product of Elliptic Curves


This package contains an implementation of the algorithm in SageMath by Ryo Yoshizumi. 

## Contents

This package provides functions to compute $(\ell,\ell)$-isogenies from a product of elliptic curves.  
The main module is `func_proposed_isogeny.py`. For the overall workflow, please refer to `test.py`.



This code is written in [SageMath](https://www.sagemath.org).



## Usage

Please select **three options** from the following:

- Method: proposed  (`"Proposed"`) or existing (`"Existing"`)
- Representation of $\ell=1^2+\cdots+1^2$ (`"One"`) or any representation (`"Square"`)
- Implementation time parameter: short (`1`) or long (`2`)  

Then, run the following command to compute and measure the runtime of an $(\ell,\ell)$-isogeny chain:

```bash
$ sage main.py {"Proposed" or "Existing"} {"One" or "Square"} {1 or 2}
 ```

Here, the base field is $\mathbb{F}_{p^2}$ and the domain is $E_0\times E_0$ where $E_0 : y^2=x^3+x$. 

The  kernel is randomly generated, and we evaluate **two** points of the form $x=(x^{(1)},0_{E_0})$ where $x^{(1)}\in E_0$.


## Example output

The output is as follows:
 ```
% sage main.py "Proposed" "One"  2
p= 11402780996313137804419565692258934141207562497476991733713707020990899136527
isogeny chain: [1321, 571, 311, 31, 31, 31, 31, 31, 31, 31, 31, 17, 17, 17, 17, 17, 17, 17, 17, 17, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 3]
ell= 1321
Time: 5856.6251418750035 sec
ell= 571
ell= 311
ell= 31
ell= 31
ell= 31
ell= 31
ell= 31
ell= 31
ell= 31
ell= 31
ell= 17
ell= 17
.
.
.
ell= 7
ell= 7
ell= 7
ell= 7
ell= 3
Total time: 6758.846125374999 sec
 ```

The above result means that it takes about 4.28 seconds to compute the $(79,79)$-isogeny and takes about 13.71 seconds to compute the whole isogeny chain. 


## Author
- **Name**: Ryo Yoshizumi 
- **Affiliation**: Kyushu University
- **Email**: yoshizumi.ryo.483@s.kyushu-u.ac.jp


