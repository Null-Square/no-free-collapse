"""Adversarial inputs, exact rational hierarchy checks, and proof boundary tests."""
from fractions import Fraction as F
from itertools import combinations
from math import comb

import numpy as np
import pytest

from no_free_collapse.hierarchy_certificate import exact_hierarchy_certificate, exact_principal_hafnians
from no_free_collapse.subhafnian_hierarchy import (
    all_principal_hafnians, averaged_step_bound, closed_form_hybrid_bound,
    iterated_hybrid_bound, normalized_step_coefficients, sharp_edge_density_upper_bound,
    spectral_diameter, subhafnian_energy,
)


@pytest.mark.parametrize('n,r', [(6,2),(10,3),(14,4)])
@pytest.mark.parametrize('family', ['balanced','complement','matching'])
def test_exact_critical_hierarchy(n,r,family):
    a=[[F(1,n) for _ in range(n)] for _ in range(n)]
    if family=='complement':
        a=[[F(i==j)-a[i][j] for j in range(n)] for i in range(n)]
    if family=='matching':
        a=[[F(1,2) if i==j else F(0) for j in range(n)] for i in range(n)]
        for i in range(0,n,2): a[i][i+1]=a[i+1][i]=F((-1)**(i//2),2)
    c=exact_hierarchy_certificate(a,r)
    assert c['equality'] and c['defect']==0 and c['spectral_interval_verified']


@pytest.mark.parametrize('n,r',[(7,2),(9,2),(11,3)])
def test_odd_ambient_dimension_rational_certificate(n,r):
    rng=np.random.default_rng(400+n)
    z=rng.integers(-2,3,size=(n,3)); q=z@z.T
    trace=int(np.trace(q))
    a=[[F(int(v),trace) for v in row] for row in q]
    assert exact_hierarchy_certificate(a,r)['defect']>=0


@pytest.mark.parametrize('n,r,aorder',[(7,2,1),(9,3,1),(11,3,2)])
def test_roos_normalized_convolution_exact_in_odd_ambient_order(n,r,aorder):
    # This is the parity-free identity needed to invoke Roos Corollary 3.1.
    rng=np.random.default_rng(700+n)
    q=rng.integers(-2,3,size=(n,n)); q=q+q.T
    h=exact_principal_hafnians(q.tolist())
    for s in combinations(range(n),2*r):
        mask=sum(1<<i for i in s)
        convolution=sum(h[tmask]*h[mask^tmask]
                        for t in combinations(s,2*aorder)
                        for tmask in [sum(1<<i for i in t)])
        assert convolution==comb(r,aorder)*h[mask]


@pytest.mark.parametrize('n',list(range(3,16)))
def test_sharp_density_including_odd_order(n):
    k=n//2
    columns=[]
    grid=np.arange(n)
    if n%2:
        if k%2: columns.append(np.ones(n)/np.sqrt(n))
        for j in range(1,k//2+1):
            columns.extend([np.sqrt(2/n)*np.cos(2*np.pi*j*grid/n),
                            np.sqrt(2/n)*np.sin(2*np.pi*j*grid/n)])
        q=np.column_stack(columns)
        p=q@q.T
    else:
        p=np.eye(n)/2
        for i in range(0,n,2): p[i,i+1]=p[i+1,i]=.5
    assert np.allclose(p@p,p,atol=1e-12)
    assert np.allclose(np.diag(p),k/n,atol=1e-12)
    off=p-np.diag(np.diag(p)); g1=np.sum(off*off)/(n*(n-1))
    theta=4*g1/spectral_diameter(p)**2
    assert np.isclose(theta,sharp_edge_density_upper_bound(n),rtol=1e-12)


@pytest.mark.parametrize('r',[True,False,1.5,-1,0])
def test_noninteger_and_out_of_domain_order_rejected(r):
    for fun in (normalized_step_coefficients,iterated_hybrid_bound,closed_form_hybrid_bound):
        with pytest.raises(ValueError): fun(np.eye(6),r)


def test_cache_is_read_only_complete_and_bound_to_matrix():
    a=np.ones((6,6))/6
    h=all_principal_hafnians(a)
    assert len(list(h.values()))==32
    with pytest.raises(TypeError): h[0]=0
    with pytest.raises(ValueError): subhafnian_energy(a,2,hafnians={0:1})
    with pytest.raises(ValueError): subhafnian_energy(2*a,2,hafnians=h)
    assert subhafnian_energy(a+np.eye(6),2,hafnians=h)==subhafnian_energy(a,2,hafnians=h)
    with pytest.raises(ValueError): all_principal_hafnians(np.zeros((21,21)))
    with pytest.raises(ValueError): all_principal_hafnians(np.empty((0,0)))


def test_exact_domain_hypotheses_and_degeneracies():
    zero=[[0]*6 for _ in range(6)]
    assert exact_hierarchy_certificate(zero,2,0,0)['equality']
    identity=np.eye(6,dtype=int).tolist()
    assert exact_hierarchy_certificate(identity,2,1,1)['equality']
    with pytest.raises(TypeError): exact_hierarchy_certificate(np.eye(6).tolist(),2)
    with pytest.raises(TypeError): exact_hierarchy_certificate([[False]*6]*6,2)
    with pytest.raises(ValueError): exact_hierarchy_certificate(identity,2,0,0)
    with pytest.raises(ValueError): exact_hierarchy_certificate(identity,True)
    with pytest.raises(ValueError): exact_hierarchy_certificate(identity,3)
    with pytest.raises(ValueError): exact_hierarchy_certificate(identity,2,1,0)
    bad=[[0]*6 for _ in range(6)]; bad[0][1]=bad[1][0]=1
    with pytest.raises(ValueError): exact_hierarchy_certificate(bad,2)
    bad[1][0]=0
    with pytest.raises(ValueError): exact_hierarchy_certificate(bad,2)
    with pytest.raises(ValueError): exact_principal_hafnians([[0]*15]*15)


def test_zero_energy_hybrid_never_divides_by_zero():
    for a in (np.zeros((10,10)),np.eye(10),np.diag(np.arange(10))):
        assert closed_form_hybrid_bound(a,3)[:2]==(0.0,0.0)
        assert averaged_step_bound(a,3)==(0.0,0.0)


def test_support_classification_exhaustively_for_unit_edge_patterns():
    # A finite independent regression for the graph step, not the universal
    # weighted equality proof: all 2^15 unweighted six-vertex support patterns.
    edges=list(combinations(range(6),2)); index={e:i for i,e in enumerate(edges)}
    edge=lambda a,b: index[tuple(sorted((a,b)))]
    products=[[(edge(a,k),edge(a,l)) for a in range(6) if a not in (k,l)]
              for k,l in edges]
    triangles=[]
    for r in combinations(range(6),3):
        rc=set(range(6))-set(r)
        triangles.append(([edge(*e) for e in combinations(r,2)],
                          [edge(*e) for e in combinations(rc,2)]))
    survivors=[]
    for mask in range(1<<15):
        if any(len({bool(mask&(1<<i)) and bool(mask&(1<<j)) for i,j in row})>1
               for row in products): continue
        if any(sum(bool(mask&(1<<i)) for i in a)!=sum(bool(mask&(1<<i)) for i in b)
               for a,b in triangles): continue
        survivors.append(mask)
    assert len(survivors)==17  # empty, complete, and 15 perfect matchings
    for mask in survivors:
        assert mask in (0,(1<<15)-1) or all(
            sum(bool(mask&(1<<edge(i,j))) for j in range(6) if j!=i)==1 for i in range(6))
