---
layout: post
title: Fast multiplicative convolution modulo a composite number
date: 2026-11-20
description: Computing the distribution of a product of residues modulo any M in O(M log M) operations, with an exact C++ implementation.
tags: math algorithms
categories: algorithms
related_posts: false
---

<!-- prettier-ignore-start -->
> ##### A note on how this was made
>
> This post, the algorithm and the code were heavily AI-generated. I have briefly checked the results, but please read it with that in mind, and let me know if you spot errors or missing accrededation to prior work. I am experimenting with ways of using AI and thought this would be cool to share :)
{: .block-warning }
<!-- prettier-ignore-end -->

Fix an integer $$M\ge2$$ and two arrays $$f$$ and $$g$$ indexed by the residues $$0,\ldots,M-1$$. Their **multiplicative convolution** is

$$
h(z)=\sum_{x,y:\,xy\equiv z\pmod M} f(x)\,g(y),\qquad z=0,\ldots,M-1.
$$

If $$f$$ and $$g$$ are the distributions of independent random residues $$X$$ and $$Y$$, then $$h$$ is the distribution of $$XY \bmod M$$. Computing $$h$$ directly takes $$M^2$$ multiplications. This post describes how to compute $$h$$ exactly in $$O(M\log M)$$ operations for **any** modulus, and benchmarks a C++ implementation.

Two special cases are well known. For a prime modulus, a discrete logarithm turns the problem into an ordinary cyclic convolution. For $$M=2^n$$, competitive-programming libraries group residues by their power of two. I could not find a public implementation for general composite $$M$$, which is the case treated here. The algebra underneath it is not new, though; see [related work](#related-work) at the end.

## Warm-up: a prime modulus

Let $$M=p$$ be prime and let $$\gamma$$ be a primitive root, so that every nonzero residue is $$\gamma^i$$ for exactly one $$i \bmod p-1$$. Since $$\gamma^i\gamma^j=\gamma^{i+j}$$, multiplying nonzero residues adds their exponents. On nonzero residues, $$h$$ is then the cyclic convolution of the sequences $$i\mapsto f(\gamma^i)$$ and $$i\mapsto g(\gamma^i)$$, which an FFT of length $$p-1$$ computes. Any product involving $$0$$ lands on $$0$$, so

$$
h(0)=\Big(\sum_x f(x)\Big)\Big(\sum_y g(y)\Big)-\sum_{z\ne0}h(z).
$$

The FFT step can also be phrased with characters, and that is the version that generalizes. A **character** of a finite abelian group $$U$$ is a homomorphism $$\chi:U\to\mathbb C^\times$$. The FFT of $$i\mapsto f(\gamma^i)$$ is, up to indexing, the list of Fourier coefficients

$$
F(\chi)=\sum_{u\in U} f(u)\,\overline{\chi(u)}
$$

over the characters of the unit group $$U=(\mathbb Z/p\mathbb Z)^\times$$. Convolution becomes pointwise multiplication: $$H(\chi)=F(\chi)G(\chi)$$.

A composite modulus brings two complications. Many residues are neither units nor zero, and a product of two nonzero residues can be zero.

## Prime powers

Let $$M=p^r$$. Every nonzero residue can be written uniquely as $$x=p^a u$$, where $$0\le a<r$$ and $$u$$ is a unit modulo $$p^{r-a}$$. (Changing the representative of $$x$$ by $$p^r$$ changes $$u$$ by $$p^{r-a}$$, so $$u$$ is only defined modulo $$p^{r-a}$$.) Call the residues with a given $$a$$ **shell** $$a$$. Writing $$U_j=(\mathbb Z/p^j\mathbb Z)^\times$$, shell $$a$$ is a copy of $$U_{r-a}$$.

For $$x=p^a u$$ and $$y=p^b v$$ we have $$xy=p^{a+b}uv$$. If $$a+b\ge r$$, the product is zero. Otherwise it lies in shell $$s=a+b$$, with unit part $$uv \bmod p^{r-s}$$.

Each $$U_j$$ is cyclic for odd $$p$$, and modulo $$2^j$$ every unit is $$\pm5^i$$. Either way, an FFT over each shell is available.

**Characters and conductors.** Reduction modulo $$p^j$$ maps $$U_r$$ onto $$U_j$$, so every character of $$U_j$$ can be read as a character of $$U_r$$. The **conductor** of a character $$\chi$$ of $$U_r$$ is the smallest $$p^{c}$$ such that $$\chi$$ comes from $$U_c$$ in this way. The trivial character has $$c=0$$. A character therefore makes sense on shell $$a$$ exactly when $$c(\chi)\le r-a$$, and on those shells we define

$$
F_a(\chi)=\sum_{u\in U_{r-a}} f(p^a u)\,\overline{\chi(u)}.
$$

Define $$G_b$$ and $$H_s$$ from $$g$$ and $$h$$ in the same way.

**The key identity.** For every shell $$s<r$$ on which $$\chi$$ is defined,

$$
H_s(\chi)=\sum_{a+b=s}F_a(\chi)\,G_b(\chi).
$$

To see this, expand $$H_s(\chi)$$ as a sum over pairs $$(p^a u,\,p^b v)$$ with $$a+b=s$$. The character is defined modulo $$p^{r-s}$$, which divides both $$p^{r-a}$$ and $$p^{r-b}$$. So $$\overline{\chi(uv)}=\overline{\chi(u)}\,\overline{\chi(v)}$$, and the sum factors.

Now collect the coefficients of each character into a polynomial in a formal variable $$t$$:

$$
F_\chi(t)=\sum_{a} F_a(\chi)\,t^a.
$$

The identity says that $$H_\chi(t)=F_\chi(t)\,G_\chi(t)$$, keeping only the degrees on which $$\chi$$ is defined, namely $$a\le r-c(\chi)$$ for nontrivial $$\chi$$ and $$a\le r-1$$ for the trivial one. For the trivial character, the discarded degrees $$\ge r$$ are exactly the mass that lands on $$0$$, which the total $$\sum_x f(x)$$ accounts for, as in the prime case. In algebraic terms, the transform is an isomorphism

$$
\mathbb C[(\mathbb Z/p^r\mathbb Z,\cdot)]\;\cong\;\mathbb C\;\times\;\mathbb C[t]/(t^{r})\;\times\prod_{\chi\ne1}\mathbb C[t]/(t^{\,r-c(\chi)+1}),
$$

where the first factor is the total mass and the second is the trivial character. Multiplication on the right is just truncated polynomial multiplication.

**Example: $$M=4$$.** The units are $$1$$ and $$3$$. The trivial character $$\chi_0$$ is identically $$1$$. The other character $$\chi_1$$ has $$\chi_1(3)=-1$$ and conductor $$4$$, so it lives only on shell $$0$$. Shell $$1$$ is the single residue $$2$$. The transform of $$f$$ is

$$
\Big(\ \textstyle\sum_x f(x),\quad \big(f(1)+f(3)\big)+f(2)\,t \bmod t^2,\quad f(1)-f(3)\ \Big).
$$

Multiplying the components gives $$\sum_z h(z)$$, then $$h(1)+h(3)$$ together with $$h(2)=\big(f(1)+f(3)\big)g(2)+f(2)\big(g(1)+g(3)\big)$$, and finally $$h(1)-h(3)$$. That is five multiplications instead of sixteen, and the four values of $$h$$ follow by solving a small linear system.

## Composite moduli

Write $$M=\prod_{i=1}^k p_i^{r_i}$$. By the Chinese remainder theorem, $$\mathbb Z/M\mathbb Z\cong\prod_i \mathbb Z/p_i^{r_i}\mathbb Z$$ as rings, so also as multiplicative monoids. The algebra above becomes a tensor product of the local ones. Concretely:

1. Reshape each array into a $$k$$-dimensional array whose $$i$$-th axis is indexed by residues modulo $$p_i^{r_i}$$.
2. Apply the local transform along each axis in turn.
3. Each block is now a polynomial in $$t_1,\ldots,t_k$$ with a separate truncation in each variable. Multiply the corresponding blocks of $$f$$ and $$g$$.
4. Invert the local transforms along each axis and reshape back.

The input does not need to factor across the CRT coordinates. Any array works.

**Example: $$M=12=4\cdot3$$.** Modulo $$3$$ the transform is three scalars: $$\sum_x f(x)$$, $$f(1)+f(2)$$ and $$f(1)-f(2)$$. Combined with the three mod-4 components, this gives $$3\times3=9$$ blocks. The three blocks that pair with the mod-4 linear polynomial have length $$2$$, and the other six are scalars. A product costs $$6\cdot1+3\cdot3=15$$ multiplications, compared with $$144$$ directly.

### Pseudocode

```text
local_forward(w, p, r):            # w is indexed by residues mod p^r
    out.total = sum(w)
    for a in 0 .. r-1:
        W = FFT over U_(r-a) of  u -> w[p^a * u]
        for each character chi of U_(r-a):
            out[chi][a] = W[chi]     # chi read as a character of U_r
    return out

local_inverse(out, p, r):
    for a in 0 .. r-1:
        W[chi] = out[chi][a] for each character chi of U_(r-a)
        w[p^a * u] = inverse FFT of W over U_(r-a)
    w[0] = out.total - (sum of all other entries of w)
    return w

multiply(f, g, M = p_1^r_1 * ... * p_k^r_k):
    reshape f and g into k-dimensional arrays by CRT
    apply local_forward along every axis of f and of g
    for each block (a choice of "total" or a character on every axis):
        multiply the two blocks as polynomials in t_1, ..., t_k,
        truncating t_i to the block's length on axis i
    apply local_inverse along every axis of the product
    reshape back to residue order
```

## Cost

**Transforms.** Along axis $$i$$, the shell FFTs touch every entry once and have lengths at most $$\varphi(p_i^{r_i})$$, so they cost $$O(M\log p_i^{r_i})$$. With [Bluestein's algorithm](https://en.wikipedia.org/wiki/Chirp_Z-transform) for awkward lengths, the total is $$O(M\log M)$$.

**Block products.** Schoolbook multiplication of truncated polynomials of length $$\ell$$ uses $$\ell(\ell+1)/2$$ products. Summing over the blocks of $$p^r$$ gives

$$
C_{p^r}=\sum_{j=1}^{r}p^j-r+1\;\le\;\frac{p^r}{1-1/p}.
$$

The counts multiply across axes, so the whole step uses at most $$M\prod_{p\mid M}(1-1/p)^{-1}=O(M\log\log M)$$ products, by Mertens' theorem. When $$M$$ is squarefree, every block is a scalar and the count is exactly $$M$$.

**The natural baseline.** A simpler method, and essentially the one used by existing $$2^n$$ implementations, transforms each shell once and then, for every pair of shells $$(a,b)$$, multiplies coefficients pointwise and adds the result into shell $$a+b$$. That includes the pairs that land on zero. I call this **shell pairs**. Locally it uses $$\sum_{j=1}^r p^j+(r^2+r+2)/2$$ products, or $$\prod_{p\mid M}(p+2)$$ for squarefree $$M$$. It is also near-linear. So the block method does not improve the asymptotics over shell pairs. It only reduces the constant in the multiplication step.

**Exact arithmetic.** The C++ implementation works over prime fields $$\mathbb F_P$$ with $$P<2^{31}$$, choosing $$P$$ so that $$P-1$$ is divisible by every transform length it needs. It bounds the largest possible $$\lvert h(z)\rvert$$ from the inputs and uses enough primes that the Chinese remainder theorem recovers the signed integer result exactly.

## Benchmarks

All times are in milliseconds, single-threaded, on an Apple M5 Max (Apple clang 21, `-O3`). For each modulus, five separate processes each timed 16 seeded dense input pairs with integer weights in $$[-7,7]$$. The table reports the median over processes of the per-process medians. Pairs are written as blocks / shell pairs.

|   $$M$$ | Factorization                                | Products, blocks / shell pairs |   Multiply step | Full transform pass | Exact, end to end | Direct |
| ------: | :------------------------------------------- | -----------------------------: | --------------: | ------------------: | ----------------: | -----: |
|   1,024 | $$2^{10}$$                                   |                  2,037 / 2,102 | 0.0017 / 0.0018 |       0.048 / 0.048 |             0.078 |   0.97 |
|   2,187 | $$3^7$$                                      |                  3,273 / 3,308 | 0.0028 / 0.0027 |       0.858 / 0.858 |             1.066 |   4.38 |
|   5,040 | $$2^4\cdot3^2\cdot5\cdot7$$                  |                10,395 / 41,328 | 0.0093 / 0.0365 |       0.301 / 0.312 |             0.361 |  23.32 |
|  11,520 | $$2^8\cdot3^2\cdot5$$                        |                27,665 / 61,264 | 0.0255 / 0.0602 |       0.829 / 0.848 |             0.901 | 122.06 |
|  30,030 | $$2\cdot3\cdot5\cdot7\cdot11\cdot13$$        |               30,030 / 245,700 | 0.0268 / 0.2071 |       2.308 / 2.310 |             2.591 |      — |
|  45,045 | $$3^2\cdot5\cdot7\cdot11\cdot13$$            |               55,055 / 196,560 | 0.0477 / 0.1677 |       2.946 / 2.949 |             3.264 |      — |
| 360,360 | $$2^3\cdot3^2\cdot5\cdot7\cdot11\cdot13$$    |            660,660 / 4,127,760 |   0.581 / 3.556 |       27.93 / 29.47 |             30.41 |      — |
| 510,510 | $$2\cdot3\cdot5\cdot7\cdot11\cdot13\cdot17$$ |            510,510 / 4,668,300 |   0.447 / 3.927 |       53.63 / 54.16 |             57.57 |      — |

The columns mean the following:

- **Multiply step** times only the loop over transformed coefficients.
- **Full transform pass** covers forward transforms, block products and inverse transforms over one prime field, with plans already built.
- **Exact, end to end** adds planning, prime selection and CRT reconstruction to the signed integer answer. One prime sufficed for every input here.
- **Direct** is the $$M^2$$ loop. I timed it only up to $$11{,}520$$.

Every trial checked that the two field methods agree bin by bin. Up to $$11{,}520$$, it also checked the exact output against direct enumeration.

Three things stand out:

- **Against direct enumeration** the gain grows with $$M$$, as expected for $$O(M\log M)$$ against $$O(M^2)$$. At $$11{,}520$$ the exact computation is about **135×** faster.
- **Against shell pairs** the block method uses up to **9× fewer products** ($$510{,}510$$ against $$4{,}668{,}300$$ at $$M=510{,}510$$), and the multiply step is about 8.8× faster. But the multiply step is a small part of the work, about 2% of the full pass for blocks and at most 12% for shell pairs. The transforms and data movement dominate, so the full pass is only 0–6% faster, with the largest gains at $$5{,}040$$ and $$360{,}360$$. For prime powers, where most characters have full conductor, the two methods are essentially the same.
- **$$2{,}187=3^7$$ is slow for its size.** Its largest shell FFT has length $$\varphi(3^7)=1{,}458$$, which the implementation handles with Bluestein's algorithm. $$5{,}040$$ needs only very short transforms.

The real gain is generality: one $$O(M\log M)$$ method, with exact integer output, for every modulus.

[Download the C++ and Python code]({{ '/assets/code/multiplicative-convolution.zip' | relative_url }}). The bundle includes a self-test against direct enumeration on 262 moduli, the benchmark driver and a NumPy prototype.

## Related work

- **The algebra.** The decomposition of $$\mathbb C[(\mathbb Z/p^r\mathbb Z,\cdot)]$$ by characters of the unit group is not new. $$\mathbb Z/p^r\mathbb Z$$ is a finite commutative chain ring, and J. A. Wood factored the semigroup determinant of such rings ("Factoring the semigroup determinant of a finite commutative chain ring", in _Coding Theory, Cryptography and Related Areas_, Springer, 2000, pp. 249–259). B. Steinberg generalized this to commutative monoids whose non-units form a principal ideal generated by a nilpotent element ([Corollary 6.13 of "Factoring the Dedekind–Frobenius determinant of a semigroup"](https://arxiv.org/pdf/2108.02166#page=32), _J. Algebra_ 605 (2022) 1–36). There, characters of the unit group and their stabilizers, which correspond to the conductors above, determine the factors. What this post adds is the algorithm built on that structure: shell FFTs, truncated block products, the CRT tensor product for composite $$M$$, and an exact implementation.
- **Prime moduli.** The discrete-log reduction is standard. It appears as the Library Checker problem [`mul_modp_convolution`](https://judge.yosupo.jp/problem/mul_modp_convolution), with implementations in [maspy's library](https://maspypy.github.io/library/mod/multiplicative_convolution_mod_p.hpp) and [cp-algorithms](https://lib.cp-algorithms.com/verify/poly/convolution_mul.test.cpp).
- **Powers of two.** [suisen's library](https://suisen-cp.github.io/cp-library-cpp/library/convolution/multiplicative_convolution_mod_2%5En.hpp) (which cites maspy and tatyam) and [maspy's library](https://maspypy.github.io/library/mod/multiplicative_convolution_mod_2n.hpp) handle $$M=2^n$$ in $$O(n2^n)$$. They use the $$\pm5^i$$ coordinates and group by pairs of 2-adic valuations, essentially the shell-pair method above.
- **Fourier transforms on semigroups.** More generally, FFTs for non-group semigroups were developed by Malandro and Rockmore for the [rook monoid](https://arxiv.org/abs/0709.4175) and by Malandro for [finite inverse semigroups](https://arxiv.org/abs/0905.1340). Note that $$(\mathbb Z/p^r\mathbb Z,\cdot)$$ is not an inverse semigroup when $$r\ge2$$.
