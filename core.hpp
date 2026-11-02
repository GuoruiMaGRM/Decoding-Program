#pragma once
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pEXFactoring.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/version.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace hhy {
using E=NTL::zz_pE;
using V=std::vector<E>;
using Mat=std::vector<V>;
using Poly=V; // coefficients in increasing q-degree; [1] represents X.
using NTL::ZZ;
inline void require(bool b,const std::string& s) { if(!b) throw std::runtime_error(s); }
inline long sz(const V& v) { return static_cast<long>(v.size()); }
inline long rows(const Mat& a) { return static_cast<long>(a.size()); }
inline E one() { return E(1); }
inline bool zero(const E& x) { return NTL::IsZero(x); }
inline void trim(V& a) { while(!a.empty()&&zero(a.back())) a.pop_back(); }
inline long degree(const Poly& a) { return sz(a)-1; }
inline ZZ integer_power(long a,long b) { ZZ z; NTL::power(z,ZZ(a),b); return z; }
inline E epow(const E& a,long b) { E z; NTL::power(z,a,b); return z; }
inline E epow(const E& a,const ZZ& b) { E z; NTL::power(z,a,b); return z; }
inline ZZ pack(const E& x,long p) {
    ZZ z(0); for(long i=NTL::deg(NTL::rep(x));i>=0;--i) z=z*p+NTL::rep(NTL::coeff(NTL::rep(x),i));
    return z;
}
inline E unpack(ZZ z,long p) {
    NTL::zz_pX a; long i=0;
    while(z>0) { long c=NTL::conv<long>(z%p); NTL::SetCoeff(a,i++,c); z/=p; }
    E x; NTL::conv(x,a); return x;
}
inline ZZ pack_modulus(const NTL::zz_pX& a,long p) {
    ZZ z(0); for(long i=NTL::deg(a);i>=0;--i) z=z*p+NTL::rep(NTL::coeff(a,i)); return z;
}
inline NTL::zz_pX unpack_modulus(ZZ z,long p) {
    NTL::zz_pX a; long i=0; while(z>0) { NTL::SetCoeff(a,i++,NTL::conv<long>(z%p)); z/=p; } return a;
}
inline std::pair<long,long> prime_power(long q) {
    require(q>=2,"q must be a prime power"); long p=2;
    while(p<=q/p&&q%p) ++p;
    if(q%p) p=q;
    long a=0,v=q; while(v%p==0) { ++a; v/=p; }
    require(v==1,"q must be a prime power"); return {p,a};
}

// Matrices use column vectors: A*x. RREF also returns pivot columns.
struct Reduction { Mat a; std::vector<long> piv; };
inline Reduction rref(Mat a,long columns=-1) {
    if(a.empty()) return {a,{}};
    long width=sz(a[0]); if(columns<0) columns=width;
    for(auto& r:a) require(sz(r)==width,"ragged matrix");
    std::vector<long> piv;
    for(long j=0;j<columns&&static_cast<long>(piv.size())<rows(a);++j) {
        long k=static_cast<long>(piv.size()),h=k;
        while(h<rows(a)&&zero(a[h][j])) ++h;
        if(h==rows(a)) continue;
        std::swap(a[k],a[h]); E inv=NTL::inv(a[k][j]);
        for(long c=j;c<width;++c) a[k][c]*=inv;
        for(long i=0;i<rows(a);++i) if(i!=k&&!zero(a[i][j])) {
            E t=a[i][j]; for(long c=j;c<width;++c) a[i][c]-=t*a[k][c];
        }
        piv.push_back(j);
    }
    return {std::move(a),std::move(piv)};
}
inline long rank(const Mat& a) { return static_cast<long>(rref(a).piv.size()); }
inline Mat transpose(const Mat& a) {
    if(a.empty()) return {}; Mat b(sz(a[0]),V(rows(a)));
    for(long i=0;i<rows(a);++i) for(long j=0;j<sz(a[i]);++j) b[j][i]=a[i][j]; return b;
}
inline V mv(const Mat& a,const V& x) {
    V z(rows(a)); for(long i=0;i<rows(a);++i) {
        require(sz(a[i])==sz(x),"matrix-vector shape");
        for(long j=0;j<sz(x);++j) z[i]+=a[i][j]*x[j];
    } return z;
}
inline Mat mm(const Mat& a,const Mat& b) {
    if(a.empty()) return {};
    require(!b.empty()&&sz(a[0])==rows(b),"matrix multiplication shape");
    Mat z(rows(a),V(sz(b[0])));
    for(long i=0;i<rows(a);++i) for(long j=0;j<rows(b);++j) if(!zero(a[i][j]))
        for(long h=0;h<sz(b[j]);++h) z[i][h]+=a[i][j]*b[j][h]; return z;
}
inline Mat inverse(const Mat& a) {
    long n=rows(a); require(n>0&&sz(a[0])==n,"inverse needs square matrix");
    Mat b(n,V(2*n));
    for(long i=0;i<n;++i) { std::copy(a[i].begin(),a[i].end(),b[i].begin()); b[i][n+i]=1; }
    auto rr=rref(b,n); require(static_cast<long>(rr.piv.size())==n,"singular matrix");
    Mat z(n,V(n)); for(long i=0;i<n;++i) std::copy(rr.a[i].begin()+n,rr.a[i].end(),z[i].begin()); return z;
}
inline std::optional<V> solve(const Mat& a,const V& b,long cols) {
    require(rows(a)==sz(b),"linear solver shape");
    if(a.empty()) return V(cols);
    Mat aug=a; for(long i=0;i<rows(a);++i) aug[i].push_back(b[i]);
    auto rr=rref(aug,cols);
    for(auto& r:rr.a) { bool nz=false; for(long j=0;j<cols;++j) nz|=!zero(r[j]); if(!nz&&!zero(r[cols])) return {}; }
    V z(cols); for(long i=0;i<static_cast<long>(rr.piv.size());++i) z[rr.piv[i]]=rr.a[i][cols]; return z;
}
inline Mat kernel(const Mat& a,long cols) {
    auto rr=rref(a); std::set<long> p(rr.piv.begin(),rr.piv.end()); Mat k;
    for(long j=0;j<cols;++j) if(!p.count(j)) {
        V v(cols); v[j]=1;
        for(long i=0;i<static_cast<long>(rr.piv.size());++i) v[rr.piv[i]]=-rr.a[i][j]; k.push_back(v);
    } return k;
}
inline Mat columns(const Mat& a,const std::vector<long>& indices) {
    Mat z(rows(a),V(indices.size())); for(long i=0;i<rows(a);++i) for(size_t j=0;j<indices.size();++j) z[i][j]=a[i][indices[j]]; return z;
}
inline V take(const V& a,const std::vector<long>& ix) { V z; for(long j:ix) z.push_back(a[j]); return z; }
inline V prefix(const V& a,long n) { require(n>=0&&n<=sz(a),"prefix length"); return V(a.begin(),a.begin()+n); }
inline V plus(V a,const V& b) { require(sz(a)==sz(b),"vector shape"); for(long i=0;i<sz(a);++i) a[i]+=b[i]; return a; }
inline V minus(V a,const V& b) { require(sz(a)==sz(b),"vector shape"); for(long i=0;i<sz(a);++i) a[i]-=b[i]; return a; }

struct Field {
    long q,p,a,M,D;
    NTL::zz_pX modulus;
    E w;
    V fq,basis;
    NTL::mat_zz_p coord_inverse;
    Field(long q_,long M_,const std::optional<ZZ>& supplied={}) :q(q_),M(M_) {
        auto pa=prime_power(q); p=pa.first; a=pa.second; D=a*M;
        require(M>0&&D<=8192,"field degree outside implementation limit");
        NTL::zz_p::init(p);
        if(supplied) { modulus=unpack_modulus(*supplied,p); require(NTL::deg(modulus)==D&&NTL::IterIrredTest(modulus),"invalid supplied modulus"); }
        else NTL::BuildIrred(modulus,D);
        NTL::zz_pE::init(modulus);
        if(a==1) w=1;
        else {
            NTL::zz_pX g; NTL::BuildIrred(g,a); NTL::zz_pEX h;
            for(long i=0;i<=a;++i) NTL::SetCoeff(h,i,E(NTL::coeff(g,i)));
            NTL::vec_zz_pE roots; NTL::FindRoots(roots,h);
            require(roots.length()==a,"subfield root count");
            w=roots[0]; for(long i=1;i<a;++i) if(pack(roots[i],p)<pack(w,p)) w=roots[i];
        }
        coordinates();
    }
    E frob(E x,long j=1) const {
        j=(j%M+M)%M; for(long i=0;i<j;++i) x=epow(x,q); return x;
    }
    void coordinates() {
        require(epow(w,q)==w,"base-field embedding not fixed by q-Frobenius");
        fq.clear(); V wb(a,one()); for(long h=1;h<a;++h) wb[h]=wb[h-1]*w;
        for(long i=0;i<q;++i) {
            long v=i; E x;
            for(long h=0;h<a;++h) { x+=E(v%p)*wb[h]; v/=p; } fq.push_back(x);
        }
        std::set<ZZ> distinct; for(auto& x:fq) distinct.insert(pack(x,p));
        require(static_cast<long>(distinct.size())==q,"base-field embedding is not injective");
        NTL::zz_pX xp; NTL::SetX(xp); E z; NTL::conv(z,xp);
        basis.assign(M,one()); for(long j=1;j<M;++j) basis[j]=basis[j-1]*z;
        NTL::mat_zz_p mat; mat.SetDims(D,D);
        for(long j=0;j<M;++j) for(long h=0;h<a;++h) {
            E x=basis[j]*wb[h]; for(long i=0;i<D;++i) mat[i][j*a+h]=NTL::coeff(NTL::rep(x),i);
        }
        NTL::zz_p determinant; NTL::inv(determinant,coord_inverse,mat);
        require(!NTL::IsZero(determinant),"F_q-coordinate basis is singular");
    }
    V coords(const E& x) const {
        if(a==1) { V out(M); for(long i=0;i<M;++i) out[i]=E(NTL::coeff(NTL::rep(x),i)); return out; }
        NTL::vec_zz_p v,z; v.SetLength(D);
        for(long i=0;i<D;++i) v[i]=NTL::coeff(NTL::rep(x),i);
        NTL::mul(z,coord_inverse,v); V out(M);
        for(long j=0;j<M;++j) { E t(1); for(long h=0;h<a;++h) { out[j]+=E(z[j*a+h])*t; t*=w; } } return out;
    }
    E from_coords(const V& v) const { require(sz(v)==M,"coordinate length"); E z; for(long i=0;i<M;++i) z+=basis[i]*v[i]; return z; }
    Mat expanded(const V& values) const { Mat z; for(auto& x:values) z.push_back(coords(x)); return z; }
    long rankq(const V& values) const { return rank(expanded(values)); }
    V independent_basis(const V& values) const {
        auto rr=rref(expanded(values)); V out;
        for(long i=0;i<static_cast<long>(rr.piv.size());++i) out.push_back(from_coords(rr.a[i])); return out;
    }
    E random(std::mt19937_64& rng) const {
        NTL::zz_pX f; std::uniform_int_distribution<long> d(0,p-1);
        for(long i=0;i<D;++i) NTL::SetCoeff(f,i,d(rng)); E x; NTL::conv(x,f); return x;
    }
    E randomq(std::mt19937_64& rng) const { return fq[std::uniform_int_distribution<long>(0,q-1)(rng)]; }
};

inline Poly padd(Poly f,const Poly& g) { f.resize(std::max(f.size(),g.size())); for(size_t i=0;i<g.size();++i) f[i]+=g[i]; trim(f); return f; }
inline Poly pscale(const E& c,Poly f) { for(auto& x:f) x*=c; trim(f); return f; }
inline Poly psub(Poly f,const Poly& g) { return padd(std::move(f),pscale(-one(),g)); }
inline Poly compose(const Poly& f,const Poly& g,const Field& K) {
    if(f.empty()||g.empty()) return {};
    Poly out(f.size()+g.size()-1); Poly twist=g;
    for(size_t i=0;i<f.size();++i) {
        if(!zero(f[i])) for(size_t j=0;j<g.size();++j) out[i+j]+=f[i]*twist[j];
        for(auto& x:twist) x=K.frob(x);
    } trim(out); return out;
}
inline E evaluate(const Poly& f,E x,const Field& K) { E z; for(auto& c:f) { z+=c*x; x=K.frob(x); } return z; }
inline V evaluations(const Poly& f,const V& x,const Field& K) { V z; for(auto& v:x) z.push_back(evaluate(f,v,K)); return z; }
inline Poly cpow(Poly f,long n,const Field& K) { Poly z{one()}; while(n) { if(n&1) z=compose(z,f,K); n/=2; if(n) f=compose(f,f,K); } return z; }
inline std::pair<Poly,Poly> right_divide(Poly f,const Poly& g,const Field& K) {
    require(!g.empty(),"zero right divisor"); Poly out(std::max(0L,sz(f)-sz(g)+1));
    while(!f.empty()&&sz(f)>=sz(g)) {
        long j=sz(f)-sz(g); E c=f.back()/K.frob(g.back(),j); out[j]+=c;
        Poly h(j+1); h[j]=c; f=psub(f,compose(h,g,K));
    } trim(out); return {out,f};
}
inline Poly exact_right(const Poly& f,const Poly& g,const Field& K) { auto qr=right_divide(f,g,K); require(qr.second.empty(),"nonzero exact right-division remainder"); return qr.first; }
inline std::optional<Poly> left_quotient(Poly f,const Poly& g,const Field& K) {
    require(!g.empty(),"zero left divisor"); Poly out(std::max(0L,sz(f)-sz(g)+1));
    while(!f.empty()&&sz(f)>=sz(g)) {
        long j=sz(f)-sz(g); E c=K.frob(f.back()/g.back(),-degree(g)); out[j]+=c;
        Poly h(j+1); h[j]=c; f=psub(f,compose(g,h,K));
    } if(!f.empty()) return {}; trim(out); return out;
}
inline Poly rgcd(Poly f,Poly g,const Field& K) { while(!g.empty()) { auto r=right_divide(f,g,K).second; f=g; g=r; } require(!f.empty(),"gcd(0,0)"); return pscale(NTL::inv(f.back()),f); }
inline Poly annihilator(const V& v,const Field& K) {
    Poly a{one()}; for(auto& x:v) { E c=evaluate(a,x,K); if(!zero(c)) a=compose(Poly{-epow(c,K.q-1),one()},a,K); } return a;
}
inline Mat moore(const V& x,long k,const Field& K) {
    Mat a(sz(x),V(k)); for(long i=0;i<sz(x);++i) { E v=x[i]; for(long j=0;j<k;++j) { a[i][j]=v; v=K.frob(v); } } return a;
}
// Ordinary polynomials are used only to construct ring elements and prime ideals.
inline Poly omul(const Poly& f,const Poly& g) {
    if(f.empty()||g.empty()) return {}; Poly z(f.size()+g.size()-1);
    for(size_t i=0;i<f.size();++i) for(size_t j=0;j<g.size();++j) z[i+j]+=f[i]*g[j]; trim(z); return z;
}
inline Poly orem(Poly f,const Poly& g) {
    require(!g.empty(),"ordinary division by zero"); while(sz(f)>=sz(g)) {
        long j=sz(f)-sz(g); E c=f.back()/g.back(); for(long i=0;i<sz(g);++i) f[i+j]-=c*g[i]; trim(f);
    } return f;
}
inline Poly opowmod(Poly f,long e,const Poly& g) {
    Poly z{one()}; while(e) { if(e&1) z=orem(omul(z,f),g); e/=2; if(e) f=orem(omul(f,f),g); } return z;
}
inline Poly ogcd(Poly f,Poly g) { while(!g.empty()) { auto r=orem(f,g); f=g; g=r; } return f.empty()?f:pscale(NTL::inv(f.back()),f); }
inline bool oirreducible(const Poly& f,const Field& K) {
    long d=degree(f); if(d<1) return false; Poly x{E(0),one()},v=x;
    for(long i=1;i<=d;++i) {
        v=opowmod(v,K.q,f);
        if(i<=d/2&&degree(ogcd(psub(v,x),f))>0) return false;
    } return psub(v,orem(x,f)).empty();
}
inline Poly oirred(long d,const Field& K,std::mt19937_64& rng) {
    for(long attempt=0;attempt<100000;++attempt) { Poly f(d+1); f[d]=1; for(long j=0;j<d;++j) f[j]=K.randomq(rng); if(oirreducible(f,K)) return f; }
    throw std::runtime_error("ordinary irreducible polynomial search exhausted");
}
inline E oevaluate(const Poly& f,const E& x) { E z; for(long i=degree(f);i>=0;--i) z=z*x+f[i]; return z; }
inline Poly monomial(long d,const E& c=one()) { Poly f(d+1); f[d]=c; trim(f); return f; }
} // namespace hhy
