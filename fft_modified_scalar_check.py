"""Full reversible-shear scalar simulation for the proposed all-exclusions basis.
Only finite scalar arithmetic; not full huge-array Fourier execution.
"""
from fft_candidate_search import Producer
from fft_new_retained_basis_audit import active_for, compile_roles
from fractions import Fraction as Q
import random


def run(h=24,seed=8918):
    c=Producer(h,4)
    ret=[(('E',j),c.results[(j,)]) for j in range(h)]
    active=active_for(c,ret)
    kk=compile_roles(c,active,ret)
    rng=random.Random(seed+h)
    v=len(c.triples)
    x={t:rng.randrange(-3,4) for t in c.triples}
    y={t:rng.randrange(-3,4) for t in c.triples}
    z=[rng.randrange(-5,6) for _ in range(kk['size'])]
    orig_z=z.copy()
    def mix(reverse=False):
        for node,ins,outs in (reversed(kk['gates']) if reverse else kk['gates']):
            pivot=ins[0]
            if reverse:
                for sl in outs[1:]:z[sl]-=z[pivot]
                if len(ins)==2:z[pivot]-=z[ins[1]]
            else:
                if len(ins)==2:z[pivot]+=z[ins[1]]
                for sl in outs[1:]:z[sl]+=z[pivot]
    def read():
        E=[z[kk['rout']['E',j]] for j in range(h)]
        a=Q(sum(E),h-3)
        out={t:a-Q(sum(E[j] for j in t),2) for t in c.triples}
        for j,(t,n,coef) in enumerate(c.pieces):out[t]+=Q(coef*z[kk['pout'][j]],2)
        return out
    mix();old=read();mix(True)
    for t,sl in kk['src'].items():z[sl]+=x[t]
    mix();new=read();mix(True)
    for t,sl in kk['src'].items():z[sl]-=x[t]
    assert z==orig_z
    assert all(new[t]-old[t]==x[t] for t in c.triples)
    print(f'PASS h={h}: {v} exact output identities, arbitrary scratch restored over {kk["size"]} roles, forward/reverse mixer checks')
if __name__=='__main__':
    run(8);run(24)
