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
p= 110564446907225951023
isogeny chain: [137, 53, 23, 23, 19, 11, 7, 7, 7]
ell= 137
Time: 11.945503083989024 sec
Currently computing the remaining isogenies.
Total time: 15.739690750022419 sec
 
 
Proposed One
the_number_of_evaluation_points: 1
p= 110564446907225951023
isogeny chain: [137, 53, 23, 23, 19, 11, 7, 7, 7]
ell= 137
Time: 11.72477579198312 sec
Currently computing the remaining isogenies.
Total time: 15.515392208006233 sec
 
 
Existing Square
the_number_of_evaluation_points: 1
p= 110564446907225951023
isogeny chain: [137, 53, 23, 23, 19, 11, 7, 7, 7]
ell= 137
Time: 17.43403587496141 sec
Currently computing the remaining isogenies.
Total time: 21.02453679102473 sec
 
 
Existing One
the_number_of_evaluation_points: 1
p= 110564446907225951023
isogeny chain: [137, 53, 23, 23, 19, 11, 7, 7, 7]
ell= 137
Time: 15.187482875015121 sec
Currently computing the remaining isogenies.
Total time: 18.834266624995507 sec
 
 
Proposed Square
the_number_of_evaluation_points: 2
p= 110564446907225951023
isogeny chain: [137, 53, 23, 23, 19, 11, 7, 7, 7]
ell= 137
Time: 13.885604707989842 sec
Currently computing the remaining isogenies.
Total time: 18.73005579202436 sec
 
 
Proposed One
the_number_of_evaluation_points: 2
p= 110564446907225951023
isogeny chain: [137, 53, 23, 23, 19, 11, 7, 7, 7]
ell= 137
Time: 13.235012207995169 sec
Currently computing the remaining isogenies.
Total time: 18.081460666959174 sec
 
 
Existing Square
the_number_of_evaluation_points: 2
p= 110564446907225951023
isogeny chain: [137, 53, 23, 23, 19, 11, 7, 7, 7]
ell= 137
Time: 24.197118834010325 sec
Currently computing the remaining isogenies.
Total time: 28.77292233298067 sec
 
 
Existing One
the_number_of_evaluation_points: 2
p= 110564446907225951023
isogeny chain: [137, 53, 23, 23, 19, 11, 7, 7, 7]
ell= 137
Time: 20.33298329199897 sec
Currently computing the remaining isogenies.
Total time: 25.161946875043213 sec

```

From the above result, for the proposed "Square" method, to compute one point, it takes about 11.95 seconds to compute the $(137,137)$-isogeny and takes about 15.84 seconds to compute the whole isogeny chain. 

## Author
- **Name**: Ryo Yoshizumi 
- **Affiliation**: Kyushu University
- **Email**: yoshizumi.ryo.483@s.kyushu-u.ac.jp


