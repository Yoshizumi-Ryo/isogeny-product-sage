# $(\ell,\ell)$-isogeny from a product of elliptic curves


This package contains the implementation of the algorithm in SageMath by Ryo Yoshizumi. 

## Content 

In this package, there is functions to compute $(\ell,\ell)$-isogeny from a product of ellptic curves. The main module is `func_proposed_isogeny.py`. For the overall flow, please refer to `test.py`.

This code is written in [SageMath](https://www.sagemath.org).



## Examples

Please choose the following three items:

- poposed method (`"Poposed"`) or existed method (`"Existed"`).
- representation of $\ell=1^2+\cdots+1^2$ (`"One"`) or any representaion (`"Square"`).
- short implementation time parameter(`1`) or long implementation time parameter(`2`)  

Then, by writing the following command, you can compute and measure runtimeone of an $(\ell,\ell)$-isogeny chain.

``` 
% sage main.py {"Proposed" or "Existed"} {"One" or "Square"} {1 or 2}
 ```

Here, the base field is $\mathbb{F}_{p^2}$ and the domain is $E_0\times E_0$ where $E_0 : y^2=x^3+x$. 

The  kernel randoml generalted, and we evaluate two points of the forms $x=(x^{(1)},0_{E_0})$ where $x^{(1)}\in E_0$.

The output is as follows:
 ```
 % sage main.py "Proposed" "One"  1 
p= 276154505650672190920223
isogeny chain: [79, 67, 29, 29, 29, 19, 11, 11, 11, 11, 7, 3, 3, 3, 3, 3, 3]
ell= 79
Time: 4.279475582996383 sec
ell= 67
ell= 29
ell= 29
ell= 29
ell= 19
ell= 11
ell= 11
ell= 11
ell= 11
ell= 7
ell= 3
ell= 3
ell= 3
ell= 3
ell= 3
ell= 3
Total time: 13.715508541004965 sec
 ```



## Author
- Name: Ryo Yoshizumi 
- Affiliation: Kyushu University
- Email: yoshizumi.ryo.483@s.kyushu-u.ac.jp


