# $(\ell,\ell)$-Isogeny from a Product of Elliptic Curves


This package contains an implementation of the algorithm in SageMath by Ryo Yoshizumi. 

## Contents

This package provides functions to compute $(\ell,\ell)$-isogenies from a product of elliptic curves.  
The main module is `func_proposed_isogeny.py` and `func_chain.py`.

This code is written in [SageMath](https://www.sagemath.org).

## Usage

Run the following command to compute and measure the runtime of an $(\ell,\ell)$-isogeny chain:

```bash
$ sage main.py
 ```

Here, the base field is $\mathbb{F}_{p^2}$ and the domain is $E_0\times E_0$ where $E_0 : y^2=x^3+x$. 

The  kernel is randomly generated, and we evaluate **one point** or **two points** of the form $x=(x^{(1)},0_{E_0})$ where $x^{(1)}\in E_0$.


## Example output

The output is as follows:
```
% sage main.py 
 
Proposed Square
the_number_of_evaluation_points: 1
p= 11402780996313137804419565692258934141207562497476991733713707020990899136527
isogeny chain: [311, 31, 31, 31, 31, 31, 31, 31, 31, 17, 17, 17, 17, 17, 17, 17, 17, 17, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 3]
ell= 311
Time: 81.4489909580152 sec
Currently computing the remaining isogenies.
Total time: 94.74064175001695 sec
 
 
Proposed One
the_number_of_evaluation_points: 1
p= 11402780996313137804419565692258934141207562497476991733713707020990899136527
isogeny chain: [311, 31, 31, 31, 31, 31, 31, 31, 31, 17, 17, 17, 17, 17, 17, 17, 17, 17, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 3]
ell= 311
Time: 76.0572325840185 sec
Currently computing the remaining isogenies.
Total time: 89.9325613330002 sec
 
 
Existing Square
the_number_of_evaluation_points: 1
p= 11402780996313137804419565692258934141207562497476991733713707020990899136527
isogeny chain: [311, 31, 31, 31, 31, 31, 31, 31, 31, 17, 17, 17, 17, 17, 17, 17, 17, 17, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 3]
ell= 311
Time: 137.05124499998055 sec
Currently computing the remaining isogenies.
Total time: 149.65957975000492 sec
 
 
Existing One
the_number_of_evaluation_points: 1
p= 11402780996313137804419565692258934141207562497476991733713707020990899136527
isogeny chain: [311, 31, 31, 31, 31, 31, 31, 31, 31, 17, 17, 17, 17, 17, 17, 17, 17, 17, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 3]
ell= 311
Time: 90.68544816601207 sec
Currently computing the remaining isogenies.
Total time: 103.47505300000194 sec
 
 
Proposed Square
the_number_of_evaluation_points: 2
p= 11402780996313137804419565692258934141207562497476991733713707020990899136527
isogeny chain: [311, 31, 31, 31, 31, 31, 31, 31, 31, 17, 17, 17, 17, 17, 17, 17, 17, 17, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 3]
ell= 311
Time: 96.02074120898033 sec
Currently computing the remaining isogenies.
Total time: 112.26134029199602 sec
 
 
Proposed One
the_number_of_evaluation_points: 2
p= 11402780996313137804419565692258934141207562497476991733713707020990899136527
isogeny chain: [311, 31, 31, 31, 31, 31, 31, 31, 31, 17, 17, 17, 17, 17, 17, 17, 17, 17, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 3]
ell= 311
Time: 83.94251741698827 sec
Currently computing the remaining isogenies.
Total time: 100.57693412501249 sec
 
 
Existing Square
the_number_of_evaluation_points: 2
p= 11402780996313137804419565692258934141207562497476991733713707020990899136527
isogeny chain: [311, 31, 31, 31, 31, 31, 31, 31, 31, 17, 17, 17, 17, 17, 17, 17, 17, 17, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 3]
ell= 311
Time: 210.08624183299253 sec
Currently computing the remaining isogenies.
Total time: 225.52677458297694 sec
 
 
Existing One
the_number_of_evaluation_points: 2
p= 11402780996313137804419565692258934141207562497476991733713707020990899136527
isogeny chain: [311, 31, 31, 31, 31, 31, 31, 31, 31, 17, 17, 17, 17, 17, 17, 17, 17, 17, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 3]
ell= 311
Time: 118.44248887500726 sec
Currently computing the remaining isogenies.
Total time: 134.0053403749771 sec

```

From the above result, for the proposed "Square" method, to compute one point, it takes about 81.45 seconds to compute the $(311,311)$-isogeny and takes about 94.74 seconds to compute the whole isogeny chain. 

## Author
- **Name**: Ryo Yoshizumi 
- **Affiliation**: Kyushu University
- **Email**: yoshizumi.ryo.483@s.kyushu-u.ac.jp


