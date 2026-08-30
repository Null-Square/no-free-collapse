# Six-variable range-only hafnian bound

This note sharpens the universal `1/192` six-variable bound by keeping track of
the asymmetric upper and lower ranges of the quadratic edge polynomial.
No PSD assumption is used.

Let `C` be a real symmetric `6 x 6` zero-diagonal matrix and write

\[
r(x)=\sum_{i<j} C_{ij}x_i x_j.
\]

Define

\[
a=-\min_x r(x),
\qquad
s=\max_x r(x).
\]

Since `r` has zero Boolean mean, `a,s>=0`.

## Theorem

For every such `C`,

\[
\boxed{
|\operatorname{haf}(C)|
\le
\frac{a\,s\,(a+s)}{48}.
}
\tag{1}
\]

When `a=s=1/2`, this reduces to

\[
|\operatorname{haf}(C)|\le\frac1{192},
\]

so the previous Chebyshev/range bound is the symmetric special case.

The point of (1) is not merely a different proof.  Combined with the exact
minimum-trace completion cost `tau(C)`, it eliminates every off-diagonal pattern
whose PSD completion penalty is not extremely small.

## Quotient by global sign

Because `r(x)=r(-x)`, quotient the six-dimensional cube by global sign.  Set

\[
y_i=x_1x_{i+1},\qquad i=1,\ldots,5.
\]

Then every pair monomial becomes either a linear or quadratic monomial in the
five independent variables `y`:

- `x_1 x_j` becomes `y_{j-1}`;
- `x_i x_j`, with `i,j>1`, becomes `y_{i-1}y_{j-1}`.

Thus there is a degree-at-most-two function `R` on `{+/-1}^5` with

\[
r(x)=R(y).
\]

The six-variable full parity becomes

\[
x_1x_2x_3x_4x_5x_6
= y_1y_2y_3y_4y_5
=:\chi(y).
\]

Since `R` has no constant term and degree at most two,

\[
\mathbb E R=0,
\qquad
\mathbb E\chi R=0.
\]

Also `R^2` has Boolean degree at most four, so

\[
\mathbb E\chi R^2=0.
\]

Let `H_+` and `H_-` be the two parity classes `chi=+1` and `chi=-1`, each
containing 16 points.  The identities above imply that the two class
populations have

\[
\mathbb E_+R=\mathbb E_-R=0
\]

and a common second moment

\[
\mathbb E_+R^2=\mathbb E_-R^2=:v.
\]

The first moment at which the two parity classes can differ is the third.
Indeed

\[
\widehat{R^3}(\{1,\ldots,5\})
=\frac12\left(\mathbb E_+R^3-\mathbb E_-R^3\right).
\]

Under the quotient this is exactly the six-variable full-parity coefficient of
`r^3`.

## Sharp third-moment interval

Let `X` be any zero-mean random variable supported in `[-a,s]` with variance

\[
\mathbb E X^2=v.
\]

For the upper third moment, the pointwise nonnegative polynomial

\[
(s-X)\left(X+\frac{v}{s}\right)^2\ge0
\]

gives, after taking expectations,

\[
\mathbb E X^3\le sv-\frac{v^2}{s}.
\]

Similarly

\[
(X+a)\left(X-\frac{v}{a}\right)^2\ge0
\]

gives

\[
\mathbb E X^3\ge -av+\frac{v^2}{a}.
\]

Therefore two zero-mean populations in the same interval with the same
variance can differ in third moment by at most

\[
(a+s)v-\left(\frac1a+\frac1s\right)v^2
=(a+s)v\left(1-\frac{v}{as}\right).
\]

The right-hand side is maximized at `v=as/2`, yielding

\[
\left|m_3^{(+)}-m_3^{(-)}\right|
\le\frac{as(a+s)}4.
\]

Hence

\[
\left|\widehat{r^3}([6])\right|
\le\frac{as(a+s)}8.
\]

Finally, the perfect-matching coefficient identity gives

\[
\widehat{r^3}([6])=3!\operatorname{haf}(C)=6\operatorname{haf}(C),
\]

and (1) follows.

## Relation to the PSD candidate

Let `tau(C)` be the minimum trace of a PSD diagonal completion, as defined in
[`diagonal_completion.md`](diagonal_completion.md).  Boolean nonnegativity
alone forces

\[
\tau(C)\ge a,
\]

because a PSD completion with trace `t` has Boolean values `t+r(x)>=0`.

The desired PSD inequality is

\[
|\operatorname{haf}(C)|
\le
\frac{\tau s(\tau+s)}{54}.
\tag{2}
\]

The range-only theorem already implies (2) whenever

\[
\frac{a s(a+s)}{48}
\le
\frac{\tau s(\tau+s)}{54},
\]

that is,

\[
\boxed{
\tau(\tau+s)\ge\frac98\,a(a+s).
}
\tag{3}
\]

Thus any counterexample to the PSD conjecture must lie in the complementary
thin shell

\[
\tau(\tau+s)<\frac98\,a(a+s).
\]

Writing `tau=a+delta`, condition (3) can fail only if

\[
\delta(2a+s+\delta)<\frac18 a(a+s).
\]

In particular,

\[
\frac{\delta}{a}
<
\frac{a+s}{8(2a+s)}
\le\frac18.
\]

So a possible counterexample must have minimum PSD completion cost less than
12.5 percent above the bare Boolean nonnegativity cost, with a stricter bound
except in the highly asymmetric limit.

## Equality geometries and the missing factor

For both known PSD equality geometries, `tau=a`:

- three disjoint pairs have `a=tau=s=1/2`;
- the equal rank-one point has `a=tau=1/6` and `s=5/6`.

At both points the actual hafnian is exactly `8/9` of the range-only upper
bound.  Therefore the remaining PSD theorem can be viewed as an `8/9`
improvement of the sharp range-only moment estimate on the boundary where PSD
adds no diagonal-completion cost.

The non-PSD `1/200` conference witness behaves differently:

\[
a=s=\frac12,
\qquad
|\operatorname{haf}(C)|=\frac1{200}
=\frac{24}{25}\cdot\frac1{192},
\]

so it lies much closer to the range-only extremum.  But its exact PSD completion
cost is

\[
\tau=\frac{3\sqrt5}{10}\approx0.67082,
\]

which is far outside the hard shell.  The generic bound (1), together with this
completion cost alone, already proves the desired PSD inequality for that
witness.

## Structural boundary `tau=a`

The shell boundary has a useful geometric interpretation.  Suppose
`tau(C)=a`, let `B` be a minimum-trace PSD completion, and let `x_*` attain
`r(x_*)=-a`.  Then

\[
x_*^TBx_*=a-a=0.
\]

Since `B` is PSD, this forces

\[
Bx_*=0.
\]

After switching coordinates by `x_*`, the minimum Boolean vertex becomes the
all-ones vector and the completion satisfies

\[
B\mathbf1=0.
\]

Thus exact zero completion gap is equivalent to a PSD completion with a Boolean
null vector.  In Gram form this means that, after a sign switch, the six Gram
vectors sum to zero.

This gives a much smaller next proof target:

> first prove the missing `8/9` improvement for minimum-trace PSD completions
> with a Boolean null vector, then control the thin perturbative shell
> `tau-a>0` using completion duality or stability.

That boundary formulation contains both known equality families.
