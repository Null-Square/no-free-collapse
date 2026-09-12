"""Independent exact-polynomial, rational, randomized, and boundary checks."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product

import numpy as np
import pytest

from no_free_collapse.zero_diagonal import (
    defect_terms, exact_certificate, recover_matching,
    recover_six_hafnian_matching, six_hafnian,
    subhafnian_energy_direct, subhafnian_energy_trace,
)


def matching(n, phase=1):
    a = np.zeros((n,n), dtype=complex if np.iscomplexobj(phase) else float)
    for i in range(0,n,2):
        a[i,i+1] = a[i+1,i] = phase
    return a


@pytest.mark.parametrize("n", range(4,13))
@pytest.mark.parametrize("complex_entries", [False, True])
def test_random_trace_and_defect_identity(n, complex_entries):
    rng = np.random.default_rng(20260912+n)
    for _ in range(50):
        a = rng.normal(size=(n,n))
        if complex_entries:
            a = a + 1j*rng.normal(size=(n,n))
        a = (a+a.T)/2
        np.fill_diagonal(a,0)
        a /= np.linalg.norm(a,2)
        terms = defect_terms(a,1)
        assert np.isclose(terms.defect, terms.term_sum, atol=2e-11)
        assert min(terms.occupancy, terms.row_variance,
                   terms.wedge_spreading, terms.spectral_slack) >= -2e-11
        assert np.isclose(subhafnian_energy_direct(a),subhafnian_energy_trace(a),atol=2e-11)


@pytest.mark.parametrize("n", [4,6,8,10,12])
@pytest.mark.parametrize("scale", [F(0),F(1,7),F(1),F(3,2)])
def test_exact_sharp_matching(n,scale):
    a = [[F(0) for _ in range(n)] for _ in range(n)]
    for i in range(0,n,2):
        a[i][i+1]=a[i+1][i]=scale*(-1 if i%4 else 1)
    result = exact_certificate(a,scale)
    assert result['defect'] == 0
    assert result['identity_residual'] == 0
    assert all(v == 0 for v in result['terms'].values())


@pytest.mark.parametrize("n", [4,5,6,8])
def test_rational_certificates(n):
    rng = np.random.default_rng(41+n)
    for _ in range(15):
        z = rng.integers(-3,4,size=(n,n));z=z+z.T;np.fill_diagonal(z,0)
        bound=int(np.max(np.sum(np.abs(z),axis=1)))
        cert=exact_certificate(z.tolist(),bound)
        assert cert['defect'] >= 0


@pytest.mark.parametrize("n", [4,6,8])
def test_phase_matching_equality(n):
    rng=np.random.default_rng(n)
    for _ in range(20):
        a=np.zeros((n,n),complex)
        perm=rng.permutation(n)
        for i in range(0,n,2):
            j,k=perm[i:i+2];a[j,k]=a[k,j]=np.exp(1j*rng.uniform(-np.pi,np.pi))
        assert abs(defect_terms(a,1).defect)<1e-11
        if n==6:
            assert abs(abs(six_hafnian(a))-1)<1e-11


@pytest.mark.parametrize("complex_entries", [False,True])
def test_stability_recovery(complex_entries):
    rng=np.random.default_rng(1234)
    for n in (6,8,10):
        m=matching(n,1j if complex_entries else -1)
        for strength in (0,1e-6,1e-5,1e-4):
            z=rng.normal(size=(n,n))
            if complex_entries:z=z+1j*rng.normal(size=(n,n))
            z=(z+z.T)/2;np.fill_diagonal(z,0)
            a=m+strength*z;a/=max(1.0,np.linalg.norm(a,2))
            rounded,upper=recover_matching(a)
            assert np.linalg.norm(a-rounded,'fro')**2 <= upper+1e-10
            assert np.allclose(rounded,m,atol=1e-2)
            if n==6:
                rounded,upper=recover_six_hafnian_matching(a)
                assert np.linalg.norm(a-rounded,'fro')**2 <= upper+1e-10


def test_half_diagonal_contraction_not_projection():
    a=np.eye(6)/2+.2*matching(6)
    assert not np.allclose(a@a,a)
    b=2*a-np.eye(6)
    terms=defect_terms(b,1)
    assert terms.subhafnian_energy/16 <= terms.edge_energy/16+1e-12


def test_rank_one_boundary_and_nonconstant_diagonal_warning():
    p=np.diag([1.,0,0,0,0,0])
    off=p-np.diag(np.diag(p))
    assert subhafnian_energy_direct(off)==0
    with pytest.raises(ValueError):defect_terms(2*p-np.eye(6),1)
    p=np.ones((6,6))/6
    b=2*(p-np.diag(np.diag(p)))
    assert np.linalg.norm(b,2)>1.6
    with pytest.raises(ValueError):defect_terms(b,1)


def test_invalid_inputs_and_singular_psd_pivots():
    for a in (np.zeros((3,3)),np.zeros((4,5)),np.eye(4),np.ones((4,4))*np.nan):
        with pytest.raises(ValueError):defect_terms(a)
    a=np.zeros((4,4));a[0,1]=1
    with pytest.raises(ValueError):defect_terms(a)
    with pytest.raises(ValueError):defect_terms(matching(4),.5)
    with pytest.raises(ValueError):defect_terms(matching(4),np.nan)
    with pytest.raises(ValueError):exact_certificate(matching(4).astype(int).tolist(),0)
    with pytest.raises(TypeError):exact_certificate(matching(4).tolist(),1)
    with pytest.raises(TypeError):exact_certificate([[False]*4]*4,1)
    with pytest.raises(ValueError):exact_certificate([[0,1,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],1)
    assert exact_certificate(matching(4).astype(int).tolist(),1)['spectral_bound_verified']
    assert exact_certificate([[0]*4 for _ in range(4)],0)['defect']==0
    with pytest.raises(ValueError):recover_matching(np.zeros((6,6)))
    with pytest.raises(ValueError):recover_six_hafnian_matching(np.zeros((6,6)))
    with pytest.raises(ValueError):six_hafnian(matching(4))
    with pytest.raises(ValueError):recover_matching(np.zeros((5,5)))


# Sparse rational polynomials: compare every monomial, not sampled points.
def add(*polys):
    out=defaultdict(F)
    for p in polys:
        for mon,c in p.items():out[mon]+=c
    return {m:c for m,c in out.items() if c}


def mul(p,q):
    out=defaultdict(F)
    for m,c in p.items():
        for k,d in q.items():out[tuple(sorted(m+k))]+=c*d
    return {m:c for m,c in out.items() if c}


def scale(p,c):return {m:v*c for m,v in p.items() if v*c}
def sq(p):return mul(p,p)
def cadd(*items):return add(*(p[0] for p in items)),add(*(p[1] for p in items))
def cmul(p,q):return add(mul(p[0],q[0]),scale(mul(p[1],q[1]),-1)),add(mul(p[0],q[1]),mul(p[1],q[0]))
def conj(p):return p[0],scale(p[1],-1)
def normsq(p):return add(sq(p[0]),sq(p[1]))


@pytest.mark.parametrize("n,complex_entries", [(4,False),(5,False),(6,False),(8,False),(4,True),(6,True)])
def test_universal_polynomial_identity(n,complex_entries):
    zero=({},{});a=[[zero for _ in range(n)] for _ in range(n)]
    count=0
    for i,j in combinations(range(n),2):
        re={(count,):F(1)};count+=1
        im={(count,):F(1)} if complex_entries else {}
        count+=int(complex_entries)
        a[i][j]=a[j][i]=(re,im)
    weights=[[normsq(x) for x in row] for row in a]
    rows=[add(*row) for row in weights]
    t=scale(add(*rows),F(1,2))
    h={}
    for i,j,k,l in combinations(range(n),4):
        term=cadd(cmul(a[i][j],a[k][l]),cmul(a[i][k],a[j][l]),cmul(a[i][l],a[j][k]))
        h=add(h,normsq(term))
    gram=[[cadd(*(cmul(a[i][k],conj(a[j][k])) for k in range(n))) for j in range(n)] for i in range(n)]
    tr4=add(*(normsq(x) for row in gram for x in row))
    trace_formula=add(scale(sq(t),F(1,2)),
                      scale(add(*(sq(x) for row in weights for x in row)),F(1,2)),
                      scale(add(*(sq(r) for r in rows)),-1),scale(tr4,F(1,4)))
    assert add(h,scale(trace_formula,-1))=={}
    ell2={(count,count):F(1)}
    occupancy=scale(mul(t,add(scale(ell2,F(n,2)),scale(t,-1))),F(n-4,2*n))
    variance=scale(add(*(sq(add(r,scale(t,F(-2,n)))) for r in rows)),F(1,2))
    wedge=add(*(mul(weights[i][j],weights[i][k]) for i in range(n)
                for j,k in combinations([v for v in range(n) if v!=i],2)))
    spectral=scale(add(scale(mul(ell2,t),2),scale(tr4,-1)),F(1,4))
    defect=add(scale(mul(ell2,t),F(n-2,4)),scale(h,-1))
    assert add(defect,scale(add(occupancy,variance,wedge,spectral),-1))=={}
