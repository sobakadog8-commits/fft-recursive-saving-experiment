"""Research-only parameter search for the attributed NStar3 finite producer.
Reconstructed from shea256/fourier-transform-below-nlogn (2026-10-08),
original producer credited to Swapnil Jain (Apache-2.0).
Modifications: reorganized program, alternate parameter checker, no upstream hashes.
Original: https://github.com/Swapnil-jain/integer-mult-kappa
Pinned commit: f2176bc1124821bf17eb63725bd366d7bdc020a3
NOT a proof of FFT bound.
"""
from collections import defaultdict, Counter
from itertools import combinations
from math import log
import sys, argparse

class Producer:
    def __init__(self,h=24,base=4):
        assert h%2==0 and h>=8
        self.h=h; self.base=base
        self.triples=list(combinations(range(h),3)); self.tid={t:i for i,t in enumerate(self.triples)}
        self.sup=[0]+[1<<i for i in range(len(self.triples))]
        self.args=[None]*len(self.sup);self.look={}; self.pieces=[];self.retained=[];self.twos={};self.top=True;self.cov={}
        w={t:self.tid[t]+1 for t in self.triples}
        self.results=self.triple(list(range(h)),w)
        for s in self.triples:self.emit(s,self.results[s],1)
        for j in range(h-1):self.retained.append((('E',j),self.results[(j,)]))
        self.retained.append((('*',),self.results[()]))
        self.star_pieces();self.prune()
    def x(self,*p):return self.tid[tuple(sorted(p))]+1
    def add(self,a,b):
        if a==0:return b
        if b==0:return a
        k=self.sup[a]|self.sup[b]
        assert not self.sup[a]&self.sup[b], 'cancellation'
        n=self.look.get(k)
        if n is None:
            n=len(self.sup);self.look[k]=n;self.sup.append(k);self.args.append((a,b))
        return n
    def total(self,vs):
        if not vs:return 0
        if len(vs)==1:return vs[0]
        k=len(vs)//2
        return self.add(self.total(vs[:k]),self.total(vs[k:]))
    def vector(self,vs):
        pre=[0]
        for x in vs:pre.append(self.add(pre[-1],x))
        suf=[0]*(len(vs)+1)
        for i in range(len(vs)-1,-1,-1):suf[i]=self.add(vs[i],suf[i+1])
        return pre[-1],[self.add(pre[i],suf[i+1]) for i in range(len(vs))]
    @staticmethod
    def groups(pts):return [pts[i:i+2] for i in range(0,len(pts),2)]
    def block(self,pts,edges,weights):
        T=self.total;A=self.add
        if len(pts)==3:
            p1,p2,p3=pts;e=lambda a,b:edges[tuple(sorted((a,b)))]
            q1=A(e(p2,p3),weights[p2]);q2=A(e(p1,p3),weights[p3]);q3=A(e(p1,p2),weights[p1])
            O={p1:A(q1,weights[p3]),p2:A(q2,weights[p1]),p3:A(q3,weights[p2])}
            tot=A(O[p1],A(q3,e(p1,p3)))
            return tot,O,{(p1,p2):weights[p3],(p1,p3):weights[p2],(p2,p3):weights[p1]}
        if len(pts)<=4:
            tot=lambda omit:T([x for p,x in edges.items() if not set(p)&set(omit)]+[x for p,x in weights.items() if p not in omit])
            return tot(()),{a:tot((a,)) for a in pts},{(a,b):tot((a,b)) for a,b in combinations(pts,2)}
        gs=self.groups(pts); ng=len(gs); e=lambda a,b:edges[tuple(sorted((a,b)))]
        coarse={(i,j):T([e(a,b) for a in gs[i] for b in gs[j]]) for i,j in combinations(range(ng),2)}
        wt={i:T([weights[a] for a in g]+[e(a,b) for a,b in combinations(g,2)]) for i,g in enumerate(gs)}
        total,outside,far=self.block(list(range(ng)),coarse,wt)
        strips={};sums={}
        for i,g in enumerate(gs):
            other=[j for j in range(ng) if j!=i]
            for a in g:
                carry=T([weights[u] for u in g if u!=a]);vals=[T([e(u,v) for u in g if u!=a for v in gs[j]]) for j in other]
                st,one=self.vector([carry]+vals)
                strips[a]={j:z for j,z in zip(other,one[1:])};sums[a]=st
        out={};single={a:A(outside[i],sums[a]) for i,g in enumerate(gs) for a in g}
        for i,g in enumerate(gs):
            for a,b in combinations(g,2):out[a,b]=outside[i]
        for i,j in combinations(range(ng),2):
            for a in gs[i]:
                left=A(far[i,j],strips[a][j])
                for b in gs[j]:
                    cross=T([e(u,v) for u in gs[i] if u!=a for v in gs[j] if v!=b])
                    out[a,b]=A(left,A(strips[b][i],cross))
        return total,single,out
    def triple(self,pts,w):
        top=self.top;self.top=False
        subs=[c for k in range(4) for c in combinations(pts,k)]
        if len(pts)<=self.base:
            return {c:self.total([v for t,v in w.items() if not set(c)&set(t)]) for c in subs}
        gs=self.groups(pts);ng=len(gs);gof={a:i for i,g in enumerate(gs) for a in g}
        cp=defaultdict(list)
        for t,v in w.items():cp[tuple(sorted({gof[a] for a in t}))].append(v)
        coarse={c:self.total(vs) for c,vs in cp.items() if not (top and len(c)==2)}
        if top:
            for i,j in [c for c in cp if len(c)==2]:
                (i0,i1),(j0,j1)=gs[i],gs[j]
                coarse[i,j]=self.add(self.add(self.x(i0,j0,i1),self.x(i0,j0,j1)),self.add(self.x(i1,j1,i0),self.x(i1,j1,j0)))
        oc=self.triple(list(range(ng)),coarse)
        one={}
        for u in pts:
            g=gof[u];others=[j for j in range(ng) if j!=g];pieces=defaultdict(list)
            for t,v in w.items():
                if u not in t:continue
                rest=[a for a in t if a!=u]
                if any(gof[a]==g for a in rest):continue
                pieces[tuple(sorted({gof[a] for a in rest}))].append(v)
            coeff={c:self.total(vs) for c,vs in pieces.items()}
            edges={c:coeff.get(c,0) for c in combinations(others,2)}
            weights={j:coeff.get((j,),0) for j in others}
            tot,single,pair=self.block(others,edges,weights)
            const=coeff.get((),0)
            one[u]={():self.add(const,tot)}
            one[u].update({(j,):self.add(const,x) for j,x in single.items()})
            one[u].update({c:self.add(const,x) for c,x in pair.items()})
        two={}
        for u,v in combinations(pts,2):
            i,j=gof[u],gof[v]
            if i==j:continue
            others=[k for k in range(ng) if k not in (i,j)]
            vals=[self.total([w.get(tuple(sorted((u,v,a))),0) for a in gs[k]]) for k in others]
            tot,leave=self.vector([w.get((u,v),0)]+vals)
            two[u,v]={():tot}|{(k,):x for k,x in zip(others,leave[1:])}
        res={}
        for ex in subs:
            E=set(ex);G=tuple(sorted({gof[a] for a in E}));surv=[u for i in G for u in gs[i] if u not in E]
            pc={():oc[G]}
            for u in surv:pc[u,]=one[u][tuple(i for i in G if i!=gof[u])]
            for u,v in combinations(surv,2):pc[u,v]=two[u,v][tuple(i for i in G if i not in (gof[u],gof[v]))]
            if len(surv)==3:
                pc[tuple(surv)]=w.get(tuple(surv),0)
                a,b,c=surv;order=[(),(a,),(b,),(a,b),(c,),(a,c),(b,c),(a,b,c)]
            else:order=[c for k in range(len(surv)+1) for c in combinations(surv,k)]
            res[ex]=self.total([pc[c] for c in order])
        return res
    def two_vec(self,u,v):
        key=(min(u,v),max(u,v))
        if key in self.twos:return self.twos[key]
        h=self.h;ng=(h+1)//2
        others=[k for k in range(ng) if k not in (u//2,v//2)]
        vals=[self.total([self.x(u,v,a) for a in range(2*k,min(2*k+2,h))]) for k in others]
        tot,leave=self.vector(vals);r={():tot}|{(k,):x for k,x in zip(others,leave)}
        self.twos[key]=r;return r
    def star_pieces(self):
        h=self.h
        for a,b in combinations(range(h),2):
            tw=self.two_vec(a,b)
            if a//2==b//2:
                for c in range(h):
                    if c//2==a//2:continue
                    S=tuple(sorted((a,b,c)))
                    self.emit(S,self.x(a,b,c^1),-1);self.emit(S,tw[(c//2,)],-1)
                continue
            pa,pb=a^1,b^1; PP=self.add(self.x(a,b,pa),self.x(a,b,pb))
            for c in range(h):
                if c in (a,b):continue
                S=tuple(sorted((a,b,c)))
                if c==pa:
                    self.emit(S,self.x(a,b,pb),-1);self.emit(S,tw[()],-1)
                elif c==pb:
                    self.emit(S,self.x(a,b,pa),-1);self.emit(S,tw[()],-1)
                else:
                    self.emit(S,self.x(a,b,c^1),-1);self.emit(S,self.add(tw[(c//2,)],PP),-1)
    def cover(self,n):
        if n not in self.cov:
            sp=self.sup[n];mask=0
            while sp:
                low=sp&-sp
                for p in self.triples[low.bit_length()-1]:mask|=1<<p
                sp^=low
            self.cov[n]=mask
        return self.cov[n]
    def emit(self,S,n,coef):
        if not n:return
        mask=sum(1<<p for p in S)
        if self.cover(n)|mask==(1<<self.h)-1:
            a,b=self.args[n];self.emit(S,a,coef);self.emit(S,b,coef)
        else:self.pieces.append((S,n,coef))
    def prune(self):
        self.active=set();st=[n for _,n,_ in self.pieces]+[n for _,n in self.retained]
        while st:
            n=st.pop()
            if n in self.active:continue
            self.active.add(n)
            if self.args[n]:st.extend(self.args[n])
        self.additions=sum(self.args[n]!=None for n in self.active)
        self.roles=self.additions+len(self.pieces)+len(self.retained)

def rank(vectors):
    d={}
    for b in vectors:
        while b:
            i=b.bit_length()-1
            if i in d:b^=d[i]
            else:d[i]=b;break
    return len(d)

def root(h,W,m):
    def psi(a):return sum(n/W*(r/m)**(1-a) for r,n in h.items() if n)
    if psi(0)>=1:return 0.,psi(0)
    l,r=0.,.1
    for _ in range(65):
        a=(l+r)/2
        if psi(a)<1:l=a
        else:r=a
    return (l+r)/2,psi(0)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--base',type=int,nargs='+',default=[3,4,5,6])
    parser.add_argument('--h',type=int,default=24)
    args=parser.parse_args()
    for base in args.base:
        try:
            c=Producer(args.h,base)
            print('CASE',dict(h=args.h,base=base,all_nodes=len(c.sup),additions=c.additions,pieces=len(c.pieces),retained=len(c.retained),roles=c.roles,unique_used=len(c.active)),flush=True)
        except Exception as e:print('FAIL',base,repr(e),flush=True)
