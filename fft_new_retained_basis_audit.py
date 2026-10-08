"""Reproducible finite audit of an alternative retained-total basis for an
attributed NStar3 producer. Adapted for research validation, not formally verified.
Requires fft_candidate_search.py in the same directory. Standard-library Python.
"""
from collections import Counter, defaultdict
from itertools import combinations
from fractions import Fraction
import sys, math, argparse, time
from fft_candidate_search import Producer,root
# Rigorous moment upper bound; all arithmetic uses Fraction.
def log_unit_upper(y, terms=18):
    assert 1 <= y <= 2
    z=(y-1)/(y+1)
    total=Fraction(0)
    power=z
    for j in range(terms):
        total+=2*power/(2*j+1)
        power*=z*z
    return total+2*power/((2*terms+1)*(1-z*z))

def log_ratio_upper(m,r):
    den=r;twos=0
    while m>=2*den:
        twos+=1;den*=2
    x=twos*log_unit_upper(Fraction(2))+log_unit_upper(Fraction(m,den))
    scale=10**14
    return Fraction(-(-x.numerator*scale//x.denominator),scale)

def exp_positive_upper(x):
    assert 0<=x<4
    # For j>=3: j! >= 6 * 4**(j-3).
    return 1+x+x*x/2+x*x*x/(6*(1-x/4))

H=24

def active_for(c, ret):
    active=set();st=[n for _,n,_ in c.pieces]+[node for _,node in ret]
    while st:
        n=st.pop()
        if n in active:continue
        active.add(n)
        if c.args[n]:st+=list(c.args[n])
    return active

def compile_roles(c,active,ret):
    users={n:[] for n in active}
    for n in sorted(active):
        if c.args[n]:
            for pos,x in enumerate(c.args[n]): users[x].append(('g',n,pos))
    for i,(_,n,_) in enumerate(c.pieces): users[n].append(('p',i))
    for name,n in ret:users[n].append(('r',name))
    edge={};src={};pout={};rout={};gates=[];size=0
    for n in sorted(active):
        if c.args[n]:
            ins=(edge[n,0],edge[n,1]); piv=ins[0]
        else:
            piv=size;size+=1;ins=(piv,);src[c.triples[n-1]]=piv
        assert len(users[n])>0,(n,'dangling')
        outs=(piv,)+tuple(range(size,size+len(users[n])-1))
        size+=len(users[n])-1
        assert len(set(ins))==len(ins) and set(ins)&set(outs)=={piv}
        gates.append((n,ins,outs))
        for u,sl in zip(users[n],outs):
            if u[0]=='g':edge[u[1],u[2]]=sl
            elif u[0]=='p':pout[u[1]]=sl
            else:rout[u[1]]=sl
    return dict(size=size,gates=gates,src=src,pout=pout,rout=rout)

class Projections:
    def __init__(self,h):
        self.h=h;self.mask=(1<<h)-1
        self.zero=tuple([0]*h);self.full=tuple(1<<i for i in range(h))
        self.rank={self.zero:0,self.full:h}
        self.complements={}
        self.diffcache={}
        self.nonnested=[]
        self.degenerate=[]
        self.nonalternating=[]
    def make(self,basis):
        basis=list(dict.fromkeys(basis))
        assert all(v.bit_count()%2==1 for v in basis),'basis must be orthonormal (odd vectors)'
        assert all((u&v).bit_count()%2==0 for u,v in combinations(basis,2)),'non-orthogonal candidate'
        cols=tuple(self.apply_basis(basis,1<<i) for i in range(self.h))
        self.rank[cols]=len(basis)
        return cols
    @staticmethod
    def apply_basis(basis,v):
        out=0
        for b in basis:
            if (b&v).bit_count()%2:out^=b
        return out
    def apply(self,P,v):
        out=0
        while v:
            b=v&-v;out^=P[b.bit_length()-1];v^=b
        return out
    def complement(self,P):
        if P not in self.complements:
            Q=tuple(p^(1<<i) for i,p in enumerate(P))
            self.rank[Q]=self.h-self.rank[P]
            self.complements[P]=Q;self.complements[Q]=P
        return self.complements[P]
    @staticmethod
    def gf2rank(arr):
        piv={}
        for a in arr:
            while a:
                d=a.bit_length()-1
                if d in piv:a^=piv[d]
                else:piv[d]=a;break
        return len(piv)
    def edge(self,A,B):
        if A==B:return 0
        key=(A,B) if id(A)<id(B) else (B,A)
        if key in self.diffcache:return self.diffcache[key]
        d=tuple(a^b for a,b in zip(A,B))
        # Ranks alone don't establish nesting: the quotient must be an actual symmetric projector.
        AinsideB=all(self.apply(B,col)==col for col in A)
        BinsideA=all(self.apply(A,col)==col for col in B)
        if not (AinsideB or BinsideA):
            self.nonnested.append((self.rank[A],self.rank[B]));raise ValueError('non-nested edge '+str((self.rank[A],self.rank[B])))
        r=self.gf2rank(d)
        assert r==abs(self.rank[A]-self.rank[B]),('difference rank mismatch',r,self.rank[A],self.rank[B])
        if r and not any(x.bit_count()%2 for x in d):
            self.degenerate.append((self.rank[A],self.rank[B]));raise ValueError('alternating residual '+str((self.rank[A],self.rank[B])))
        self.diffcache[key]=r
        return r

def frames(c,active):
    p=Projections(c.h);triples=[sum(1<<i for i in t) for t in c.triples]
    ls=[p.make([t]) for t in triples]
    labels={}
    for n in sorted(active):
        sp=c.sup[n];common=p.mask;cover=0;ids=[]
        while sp:
            b=sp&-sp;j=b.bit_length()-1;t=triples[j]
            ids.append(j);common&=t;cover|=t;sp^=b
        if len(ids)==1 or common.bit_count()>=2:
            basis=[triples[j] for j in ids]
        else:basis=[1<<j for j in range(c.h) if cover>>j&1]
        labels[n]=p.make(basis)
    return p,labels,ls

def hist(c,kk,active,ret,verbose=False):
    h=c.h;v=len(c.triples);R=kk['size']
    p,labels,lines=frames(c,active)
    pieces={t:[] for t in c.triples}
    for j,sl in kk['pout'].items():pieces[c.pieces[j][0]].append(sl)
    retained=set(kk['rout'].values())
    assert len(retained)==h
    ops=[]
    def mix(mode,reverse=False):
        gates=reversed(kk['gates']) if reverse else kk['gates']
        ops.append([(tuple(set(ins+outs)),p.zero if mode=='low' else p.full if mode=='high' else labels[n]) for n,ins,outs in gates])
    def scat():ops.append([((sl,),p.zero) for sl in sorted(retained)])
    def inj(inv):ops.append([(tuple(pieces[t]),p.complement(lines[j]) if inv else p.zero) for j,t in enumerate(c.triples)])
    def src(high):ops.append([((kk['src'][t],),p.full if high else lines[j]) for j,t in enumerate(c.triples)])
    mix('low');scat();inj(False);mix('low',True);src(False);mix('label');scat();inj(True);mix('high',True);src(True)
    reverse_ops=[[(roles,p.complement(frame)) for roles,frame in reversed(op)] for op in reversed(ops)]
    local=[];checks=0
    for stage,actions,designated in ((1,ops,6),(2,reverse_ops,3)):
        cur=[p.zero]*R;count=Counter();pending={}
        for oi,op in enumerate(actions):
            for roles,frame in op:
                for role in roles:
                    if role in retained and oi==designated:
                        if stage==1:
                            r=p.edge(cur[role],frame)
                            if r:count[r]+=1
                        else:pending[role]=frame
                        checks+=1;continue
                    if role in pending:
                        r=p.edge(frame,pending.pop(role))
                        if r:count[r]+=1
                        checks+=1
                    r=p.edge(cur[role],frame)
                    if r:count[r]+=1
                    cur[role]=frame;checks+=1
        assert not pending
        for old in cur:
            r=p.edge(old,p.full)
            if r:count[r]+=1
            checks+=1
        local.append(count)
    assert local[0]==local[1],[(local[0]-local[1]),(local[1]-local[0])]
    ell=sum(p.rank[labels[node]] for _,node in ret)
    assert sum(r*n for r,n in local[0].items())==h*R+ell, (sum(r*n for r,n in local[0].items()),h*R+ell)
    H=Counter({r:2*v*n for r,n in local[0].items()})
    H[h*h-h]+=2*v*R
    H[(h-1)**2]+=2*v*v
    H[h-1]+=4*v*v
    H[1]+=v*v
    W=2*v*v+2*v*R
    a,phi0=root(H,W,h*h)
    out=dict(R=R,W=W,ell=ell,H=H,s=sum(k*n for k,n in H.items()),gap=h*h*W-sum(k*n for k,n in H.items()),root=a,phi0=phi0,checks=checks)
    if verbose:print('PROJECTOR_EDGES',checks,'distinct cache',len(p.diffcache),flush=True)
    return out

def coeff_check(c,alt):
    # Entire readout identity tested exhaustively, no floating-point operations
    v=len(c.triples)
    rows={t:[0]*v for t in c.triples}
    for s,node,coef in c.pieces:
        sp=c.sup[node]
        while sp:
            low=sp&-sp;rows[s][low.bit_length()-1]+=coef;sp^=low
    for s in c.triples:
        st=set(s)
        for j,t in enumerate(c.triples):
            if alt:
                # twice the readout: 2/21 *sum E_j - sum_{j in S} E_j
                # sum E_j =21 for every source triple, hence simplifies identically
                ret=2-len(st-set(t))
            else:ret= len(st&set(t))-1
            assert rows[s][j]+ret==2*(s==t),(s,t,rows[s][j],ret)
    return v*v

if __name__=='__main__':
    pp=argparse.ArgumentParser()
    pp.add_argument('--h',type=int,default=24)
    pp.add_argument('--mode',choices=['both','baseline','alternative'],default='both')
    flags=pp.parse_args()
    c=Producer(flags.h,4)
    baseline=c.retained
    alternate=[(('E',j),c.results[(j,)]) for j in range(flags.h)]
    for name,ret in [('baseline',baseline),('alternative',alternate)]:
        if flags.mode not in ('both',name): continue
        a=active_for(c,ret)
        kk=compile_roles(c,a,ret)
        print('TEST',name,'additions',sum(c.args[n]!=None for n in a),'roles',kk['size'],flush=True)
        coeff_check(c,name=='alternative')
        result=hist(c,kk,a,ret,verbose=True)
        print('RESULT',name,{k:v for k,v in result.items() if k!='H'},flush=True)
        if name=='baseline':
            expected={1:184811440,2:107199136,3:50583808,4:53478128,5:33031680,6:27611408,7:14767104,8:20195472,9:6217728,10:14119424,11:5051904,12:13742960,13:4663296,14:14864256,15:4849504,16:13014320,17:10196912,18:33371712,19:13601280,20:22774048,22:44038192,23:16479408,24:4048,529:8193152,552:199193984}
            if flags.h==24: assert dict(result['H'])==expected, (result['H']-Counter(expected),Counter(expected)-result['H'])
            if flags.h==24:print('BASELINE HISTOGRAM REPRODUCED EXACTLY',flush=True)
        else:
            # certified_psi_upper is configured for original W; adapt for new W
            def certified_upper(a):
                m=flags.h*flags.h
                return sum(Fraction(r*n,result['W']*m)*exp_positive_upper(a*log_ratio_upper(m,r)) for r,n in result['H'].items() if n)
            for a in (Fraction(74,10**6),Fraction(7402,10**8),Fraction(74025,10**9),Fraction(74026,10**9),Fraction(75,10**6)):
                bound=certified_upper(a)
                print('CERTIFICATE',str(a),float(1-bound),bound<1,flush=True)
            print('ALT_HIST',dict(sorted(result['H'].items())),flush=True)
