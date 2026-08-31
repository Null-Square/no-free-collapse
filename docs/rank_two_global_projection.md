# Global rank-two projection gradient theorem

This note closes the six-variable hafnian-gradient contraction for every real rank-two orthogonal projection.

For a real `6 x 6` orthogonal projection `P`, define

\[
q_1(P)=\sum_{i<j}P_{ij}^2,
\qquad
q_2(P)=\sum_{i<j}\operatorname{haf}(P_{\widehat i,\widehat j})^2.
\]

The theorem is

\[
\boxed{q_2(P)\le \frac14 q_1(P)}
\tag{1}
\]

for every rank-two `P`.

The proof combines three previously established pieces:

1. the balanced-diagonal theorem when every `P_ii<=1/2`;
2. the quadratic-dual high-pair theorem when the two largest diagonal entries sum to at least `3/2`;
3. the middle-strip argument proved below.

Rank four then follows immediately from projection complementation because `offdiag(I-P)=-offdiag(P)`, so `q1(I-P)=q1(P)` and every complementary `4 x 4` hafnian is unchanged (hafnian degree two).

## 1. Diagonal-only quadratic-dual certificate

Write

\[
d_i=P_{ii},\qquad 0\le d_i\le1,\qquad \sum_i d_i=2,
\]

and define

\[
r_i=d_i(1-d_i),\qquad R=\sum_i r_i,\qquad E=\sum_i r_i^2.
\]

The exact Pluecker defect plus chord linearization and the quadratic node potential

\[
\alpha(t)=-t^2+\frac32t-\frac14
\]

give the already-proved lower bound

\[
8\left(\frac14q_1-q_2\right)\ge H(d),
\tag{2}
\]

where

\[
H(d)=H_0(d)-C_{\rm low}(d)-C_{\rm high}(d),
\qquad
H_0=R(1-R)+2E,
\tag{3}
\]

with

\[
C_{\rm low}
=\sum_{d_i+d_j<1/2}
4d_i d_j(1-d_i-d_j)(1-2d_i-2d_j),
\tag{4}
\]

and

\[
C_{\rm high}
=\sum_{d_i+d_j>3/2}
4(1-d_i)(1-d_j)(d_i+d_j-1)(2d_i+2d_j-3).
\tag{5}
\]

The high-pair theorem already proves `H>=0` if some pair has sum at least `3/2`.

Thus sort

\[
d_1\ge d_2\ge\cdots\ge d_6
\]

and consider only the middle strip

\[
\boxed{d_1\ge\frac12,\qquad d_1+d_2\le\frac32.}
\tag{6}
\]

There are no high corrections in this strip, so `H=H0-C_low`.

## 2. Uniform lower bounds for the base term

Since `sum d_i=2`,

\[
R=2-\sum_i d_i^2.
\]

Cauchy gives

\[
E\ge\frac{R^2}{6},
\]

and therefore

\[
H_0\ge R-\frac23R^2.
\tag{7}
\]

We next locate `R` in the middle strip.

For the upper bound on `R`, Cauchy on the five coordinates below `d1` gives

\[
\sum_i d_i^2
\ge d_1^2+\frac{(2-d_1)^2}{5}
\ge\frac7{10},
\]

because `d1>=1/2`. Hence

\[
R\le\frac{13}{10}.
\tag{8}
\]

For the lower bound on `R`, every `d_i` for `i>=2` is at most `d2`, so

\[
\sum_i d_i^2
\le d_1^2+d_2(2-d_1).
\]

Using `d2<=min(d1,3/2-d1)` shows the right-hand side is at most `3/2`; therefore

\[
R\ge\frac12.
\tag{9}
\]

The function `R-(2/3)R^2` is concave, so on `[1/2,13/10]` its minimum occurs at an endpoint. Equations (7)--(9) yield

\[
\boxed{H_0\ge\frac{13}{75}.}
\tag{10}
\]

There is a stronger bound when `d1>=3/4`. The same Cauchy estimate gives

\[
\sum_i d_i^2\ge
\left(\frac34\right)^2+
\frac{(5/4)^2}{5}=\frac78,
\]

so `R<=9/8`; together with `R>=1/2`,

\[
\boxed{d_1\ge\frac34\quad\Longrightarrow\quad H_0\ge\frac9{32}.}
\tag{11}
\]

## 3. A universal bound for one low edge

For a low pair `x+y=s<1/2`, define

\[
c(x,y)=4xy(1-s)(1-2s).
\]

Since `4xy<=s^2`,

\[
c(x,y)\le f(s):=s^2(1-s)(1-2s).
\]

The nonzero critical point of `f` in `[0,1/2]` is

\[
s_*=\frac{9-\sqrt{17}}{16},
\]

and

\[
f(s_*)=\frac{-107+51\sqrt{17}}{4096}
<\boxed{\frac{13}{512}}.
\tag{12}
\]

Thus every low edge costs less than `13/512`.

If there are at most six low edges, then by (10) and (12),

\[
H\ge\frac{13}{75}-6\frac{13}{512}
=\frac{403}{19200}>0.
\tag{13}
\]

If `d1>=3/4`, coordinate `1` cannot belong to a low edge, so there are at most ten low edges. Equations (11)--(12) give

\[
H\ge\frac9{32}-10\frac{13}{512}
=\boxed{\frac7{256}}>0.
\tag{14}
\]

Hence the only unresolved possibility inside (6) has

\[
\frac12\le d_1<\frac34
\tag{15}
\]

and at least seven low edges among the remaining five coordinates.

## 4. Eight low edges are impossible

Write the remaining five coordinates as

\[
y_1\ge y_2\ge y_3\ge y_4\ge y_5,
\qquad
S=\sum_{i=1}^5y_i=2-d_1>\frac54.
\tag{16}
\]

Low edges form a threshold/Ferrers graph because if `y_i+y_j<1/2`, every pair farther down the sorted order is also low.

If at least eight of the ten pairs were low, at most two pairs would be non-low. Therefore

\[
y_1+y_4<\frac12
\]

(otherwise `(1,2),(1,3),(1,4)` would all be non-low), and

\[
y_2+y_3<\frac12
\]

(otherwise `(1,2),(1,3),(2,3)` would all be non-low).

Also `(y_4,y_5)` is low, so `y_5<1/4`. Consequently

\[
S=(y_1+y_4)+(y_2+y_3)+y_5<1+\frac14=\frac54,
\]

contradicting (16).

Thus exactly seven low edges remain.

## 5. The seven-edge graph is unique

Seven low edges mean exactly three non-low edges among the five sorted coordinates.

A threshold graph with exactly three non-low edges has only two possible shapes:

- the triangle `(1,2),(1,3),(2,3)`;
- the star `(1,2),(1,3),(1,4)`.

The star is impossible under `S>5/4`: it forces

\[
y_2+y_3<\frac12,\qquad y_1+y_5<\frac12,
\]

and then `y4<1/4`, giving again `S<5/4`.

Therefore the non-low edges are exactly the triangle on the three largest coordinates. Write

\[
a\ge b\ge c\ge p\ge q.
\]

The seven low edges are the six cross edges between `{a,b,c}` and `{p,q}`, together with `(p,q)`.

## 6. Exact maximum of the seven-edge correction

For one low edge let

\[
c(x,z)=4xz(1-x-z)(1-2x-2z).
\]

On the low triangle `x+z<=1/2`,

\[
\frac{\partial^2c}{\partial x^2}
=8z(6x+4z-3)\le0,
\]

and symmetrically `c` is concave in `z` for fixed `x`.

Fix the top and bottom masses

\[
X=a+b+c,\qquad Z=p+q.
\]

Jensen first equalizes the top coordinates and then the bottom coordinates. The bottom-bottom term is also maximized at `p=q` because, for fixed `Z`, it is proportional to `pq`. Hence

\[
C_{\rm low}\le C_7(S,Z),
\qquad S=X+Z,
\]

where

\[
C_7(S,Z)
=6c\left(\frac{S-Z}{3},\frac Z2\right)
+c\left(\frac Z2,\frac Z2\right).
\tag{17}
\]

The exact polynomial is

\[
C_7(S,Z)
=\frac Z9\left(
8S^3-36S^2-6SZ^2+18SZ+36S
+16Z^3-9Z^2-27Z
\right).
\tag{18}
\]

Its `S` derivative is

\[
\frac{\partial C_7}{\partial S}
=\frac{2Z}{3}
\left(4S^2-12S+6-Z^2+3Z\right)<0
\tag{19}
\]

for `1<=S<=3/2` and `0<=Z<=1/2`.

The top-top pairs must be non-low, so `(S-Z)/3>=1/4`, equivalently

\[
S\ge Z+\frac34.
\]

We may enlarge the actual domain to `S>=1`. Since (19) is negative, the maximizing boundary is

\[
S_{\min}=\max\left(1,Z+\frac34\right).
\]

For `0<=Z<=1/4`, substitute `S=1` into (18):

\[
C_7(1,Z)=\frac{Z(Z-1)(16Z^2+Z-8)}9.
\]

Its derivative is positive on `[0,1/4]` because the derivative numerator

\[
64Z^3-45Z^2-18Z+8
\]

is decreasing there and equals `27/16>0` at `Z=1/4`.

For `1/4<=Z<=1/2`, substitute `S=Z+3/4`:

\[
C_7\left(Z+\frac34,Z\right)
=\frac{Z(2Z-1)(8Z^2-2Z-9)}8.
\]

Its derivative is negative on `[1/4,1/2]` because the derivative numerator

\[
64Z^3-36Z^2-32Z+9
\]

is decreasing there and already equals `-1/4` at `Z=1/4`.

Therefore the unique maximum occurs at

\[
S=1,\qquad Z=\frac14,
\]

and

\[
\boxed{C_{\rm low}\le\frac9{64}.}
\tag{20}
\]

Combining (10) and (20),

\[
H\ge\frac{13}{75}-\frac9{64}
=\boxed{\frac{157}{4800}}>0.
\tag{21}
\]

This closes the final seven-edge case and hence the entire middle strip.

## 7. Global rank-two and rank-four conclusions

The diagonal of a rank-two projection falls into one of three regimes:

1. `d1<=1/2`: balanced-diagonal theorem;
2. `d1+d2>=3/2`: high-pair quadratic-dual theorem;
3. `d1>=1/2` and `d1+d2<=3/2`: middle-strip theorem above.

Thus (1) holds for every rank-two projection.

For `Q=I-P`,

\[
q_1(Q)=q_1(P),\qquad q_2(Q)=q_2(P),
\]

so the same theorem holds for rank four.

Together with the already-proved ranks one, three and five (and trivial ranks zero and six), we obtain the full projection theorem:

\[
\boxed{
q_2(P)\le\frac14q_1(P)
\quad\text{for every real }6\times6\text{ orthogonal projection }P.
}
\tag{22}
\]

This proves the hafnian-gradient contraction on the entire projection locus. It does **not** yet prove the analogous inequality for an arbitrary PSD contraction `0<=A<=I`; that extension is the next separate problem.
