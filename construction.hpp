#pragma once
#include "core.hpp"

namespace hhy {
struct Parameters {
    long q=2,M=16,R=4,r=2,b=1,ell=3,s=1;
    long n() const { return ell*R; }
    long k() const { return s*r+b; }
    long kg() const { return s*R+b; }
    long tau() const { return (n()-kg())/2; }
    long u() const { return (R-r)/2; }
    void check() const {
        require(1<=b&&b<=r&&r<R&&0<=s&&s<ell&&n()<=M,"invalid code parameters");
        prime_power(q);
    }
};
struct Code {
    Parameters p;
    std::string family,label;
    std::unique_ptr<Field> field;
    V beta;
    Mat gamma,polys,G;
    std::vector<Mat> local; // local[g]: r by k, maps message to local coefficients.
    Mat info_inverse;
    std::vector<long> info_indices;
    std::map<std::string,std::string> audit;
    V encode(const V& msg) const { return mv(transpose(G),msg); }
    bool member(const V& word) const {
        if(sz(word)!=p.n()) return false;
        return encode(mv(info_inverse,take(word,info_indices)))==word;
    }
    bool certify(const V& y,const V& c) const {
        return sz(c)==p.n()&&member(c)&&field->rankq(minus(y,c))<=p.tau();
    }
    void finish() {
        p.check(); Field& K=*field;
        require(K.q==p.q&&K.M==p.M,"field and code parameters differ");
        require(sz(beta)==p.n()&&K.rankq(beta)==p.n(),"evaluation points are dependent");
        require(rows(polys)==p.k()&&rows(gamma)==p.ell,"message/gamma shape");
        G.clear();
        for(auto& f:polys) { require(sz(f)<=p.kg(),"message degree exceeds ambient bound"); G.push_back(evaluations(f,beta,K)); }
        require(rank(G)==p.k(),"generator rank mismatch");
        info_indices=rref(G).piv;
        info_inverse=inverse(transpose(columns(G,info_indices)));
        local.clear();
        for(long g=0;g<p.ell;++g) {
            require(sz(gamma[g])==p.R&&K.rankq(gamma[g])==p.R,"local points dependent");
            Mat inv=inverse(moore(prefix(gamma[g],p.r),p.r,K)); Mat C;
            for(auto& row:G) {
                V part(row.begin()+g*p.R,row.begin()+(g+1)*p.R);
                V f=mv(inv,prefix(part,p.r));
                require(evaluations(f,gamma[g],K)==part,"local Gabidulin inclusion failed"); C.push_back(f);
            }
            local.push_back(transpose(C));
        }
    }
};
inline void write_vec(std::ostream& out,const V& v,long p) {
    out<<v.size(); for(auto& x:v) out<<' '<<pack(x,p); out<<'\n';
}
inline V read_vec(std::istream& in,long p) {
    long n; in>>n; require(in.good()&&n>=0&&n<=100000,"invalid vector record"); V v(n);
    for(auto& x:v) { ZZ z; in>>z; require(in.good()&&z>=0,"invalid field element"); x=unpack(z,p); } return v;
}
inline void save_code(const Code& C,const std::filesystem::path& path) {
    std::ofstream out(path); require(out.good(),"cannot write code snapshot"); auto p=C.p;
    out<<"HHYCODE 2\n"<<C.family<<' '<<C.label<<'\n';
    out<<p.q<<' '<<p.M<<' '<<p.R<<' '<<p.r<<' '<<p.b<<' '<<p.ell<<' '<<p.s<<'\n';
    out<<pack_modulus(C.field->modulus,C.field->p)<<' '<<pack(C.field->w,C.field->p)<<'\n';
    write_vec(out,C.beta,C.field->p);
    for(auto& v:C.gamma) write_vec(out,v,C.field->p);
    for(auto& v:C.polys) write_vec(out,v,C.field->p);
    for(auto& v:C.G) write_vec(out,v,C.field->p);
    out<<C.audit.size()<<'\n'; for(auto& kv:C.audit) out<<kv.first<<' '<<std::quoted(kv.second)<<'\n';
    require(out.good(),"failed to write code snapshot");
}
inline Code load_code(const std::filesystem::path& path) {
    std::ifstream in(path); require(in.good(),"cannot read code snapshot: "+path.string());
    std::string magic; long version; in>>magic>>version; require(magic=="HHYCODE"&&(version==1||version==2),"unknown code format");
    Code C; in>>C.family>>C.label;
    auto& p=C.p; in>>p.q>>p.M>>p.R>>p.r>>p.b>>p.ell>>p.s; p.check();
    ZZ modulus,w; in>>modulus>>w; C.field=std::make_unique<Field>(p.q,p.M,modulus);
    C.field->w=unpack(w,C.field->p); C.field->coordinates();
    C.beta=read_vec(in,C.field->p);
    for(long g=0;g<p.ell;++g) C.gamma.push_back(read_vec(in,C.field->p));
    for(long j=0;j<p.k();++j) C.polys.push_back(read_vec(in,C.field->p));
    Mat expected;
    if(version==2) for(long j=0;j<p.k();++j) expected.push_back(read_vec(in,C.field->p));
    long count; in>>count; require(in.good()&&count>=0&&count<1000,"bad audit record");
    for(long j=0;j<count;++j) { std::string k,v; in>>k>>std::quoted(v); C.audit[k]=v; }
    C.finish();
    if(version==2) require(C.G==expected,"snapshot generator differs from polynomial evaluation");
    return C;
}
inline Mat qkernel(const Poly& f,const Field& K) {
    Mat output=transpose(K.expanded(evaluations(f,K.basis,K)));
    return kernel(output,K.M);
}
inline V root_basis(const Poly& f,const Field& K) {
    V out; for(auto& v:qkernel(f,K)) out.push_back(K.from_coords(v));
    for(auto& x:out) require(zero(evaluate(f,x,K)),"torsion root check failed"); return out;
}
inline Mat message_cut(const std::vector<Poly>& Q,const Parameters& p,const Field& K) {
    long width=p.s*p.R+p.r; Mat A;
    for(auto& f:Q) for(long j=0;j<p.r;++j) {
        Poly v=compose(monomial(j),f,K); v.resize(width); A.push_back(v);
    }
    require(rank(A)==(p.s+1)*p.r,"HHY operator independence failed");
    Mat constraints;
    for(long j=p.kg();j<width;++j) { V row; for(auto& v:A) row.push_back(v[j]); constraints.push_back(row); }
    Mat combinations=kernel(constraints,rows(A));
    Mat result=mm(combinations,A);
    auto rr=rref(result); require(static_cast<long>(rr.piv.size())==p.k(),"message cut dimension failed");
    result.clear();
    for(long i=0;i<p.k();++i) {
        for(long j=p.kg();j<width;++j) require(zero(rr.a[i][j]),"message cut tail is nonzero");
        result.push_back(prefix(rr.a[i],p.kg()));
    } return result;
}

// HHY degree-two coefficient ring: u=1/rho, v=T/rho,
// w=T^2/rho=1-rho[1]*v-rho[0]*u. Multiplication is composition.
inline Code build_hhy(Parameters target,long source_M,long max_M,long period_cap,
                      uint64_t construction_seed,const std::string& label,
                      bool strict_M=true,long theta_power=1) {
    target.check(); require(target.ell<=target.q+1,"HHY rational-point limit ell<=q+1");
    long powerR=target.R; while(powerR%target.q==0) powerR/=target.q;
    require(powerR==1&&target.R>=target.q,"automatic HHY family requires R=q^mu, mu>=1");
    require(source_M>=2&&source_M%2==0,"HHY source degree must be even");
    std::cerr<<label<<": source module\n";
    Field B(target.q,source_M); std::mt19937_64 rng(construction_seed);
    Poly rho=oirred(2,B,rng);
    NTL::zz_pEX rhoE; for(long j=0;j<=2;++j) NTL::SetCoeff(rhoE,j,rho[j]);
    NTL::vec_zz_pE rho_roots; NTL::FindRoots(rho_roots,rhoE);
    require(rho_roots.length()==2,"HHY deleted-point roots"); E eta=rho_roots[0];
    if(pack(rho_roots[1],B.p)<pack(eta,B.p)) eta=rho_roots[1];
    NTL::zz_pX xp; NTL::SetX(xp); E theta0; NTL::conv(theta0,xp);
    theta0=epow(theta0,theta_power);
    require(B.frob(theta0)!=theta0,"SOURCE_SELECTOR: theta lies in F_q");
    E theta; bool found=false;
    for(auto& shift:B.fq) if(!zero(oevaluate(rho,theta0+shift))) { theta=theta0+shift; found=true; break; }
    require(found&&B.frob(theta)!=theta,"SOURCE_SELECTOR: no legal HHY characteristic parameter; use source_M=4");
    E aa=NTL::inv(theta-eta),cc=NTL::inv(theta-B.frob(eta));
    Poly U{aa*cc,cc+epow(aa,B.q-1)*cc,cc/aa};
    Poly Vop{theta*aa*cc,theta*cc+B.frob(eta)*epow(aa,B.q-1)*cc,B.frob(eta)*cc/aa};
    U=cpow(U,target.R,B); Vop=cpow(Vop,target.R,B);
    E a=rho[1],c=rho[0];
    require(psub(compose(U,Vop,B),compose(Vop,U,B)).empty(),"HHY source operators do not commute");
    Poly relation=psub(padd(padd(compose(Vop,Vop,B),pscale(a,compose(U,Vop,B))),pscale(c,compose(U,U,B))),U);
    require(relation.empty(),"HHY source relation failed");
    E characteristic=Vop[0]/U[0];
    require(B.frob(characteristic)!=characteristic&&!zero(oevaluate(rho,characteristic)),"invalid characteristic");
    Poly W=psub(psub(Poly{one()},pscale(a,Vop)),pscale(c,U));
    auto image=[&](const Poly& numerator,long h,const Poly& u,const Poly& v,const Poly& w) {
        require(degree(numerator)<=2*h,"ring numerator degree"); Poly out;
        for(long i=0;i<sz(numerator);++i) if(!zero(numerator[i])) {
            long j=i/2,parity=i%2; require(h-j-parity>=0,"negative ring exponent");
            auto term=compose(compose(cpow(w,j,B),cpow(v,parity,B),B),cpow(u,h-j-parity,B),B);
            out=padd(out,pscale(numerator[i],term));
        } return out;
    };
    auto phi=[&](const Poly& num,long h) { return image(num,h,U,Vop,W); };
    long nu=(std::max(2L,source_M)+target.s)/2+1,e=2*nu-target.s;
    Poly pi=oirred(e,B,rng); long h=(e+1)/2,v=2*h-e;
    Poly PI=rgcd(phi(pi,h),phi(omul(monomial(v),pi),h),B);
    require(degree(PI)==target.R*e&&!zero(PI[0]),"auxiliary ideal-isogeny degree");
    Poly UI=exact_right(compose(PI,U,B),PI,B),VI=exact_right(compose(PI,Vop,B),PI,B);
    Poly WI=psub(psub(Poly{one()},pscale(a,VI)),pscale(c,UI));
    std::vector<Poly> Q,torsion;
    for(long i=0;i<=target.s;++i) {
        Q.push_back(exact_right(phi(omul(pi,monomial(i)),nu),PI,B));
        require(degree(Q.back())==target.s*target.R,"HHY Q degree");
    }
    for(long g=0;g<target.ell;++g) {
        Poly f;
        if(g==0) f=rgcd(UI,VI,B);
        else { E z=B.fq[g-1]; f=rgcd(psub(VI,pscale(z,UI)),psub(WI,pscale(z,VI)),B); }
        require(degree(f)==target.R&&!zero(f[0]),"rational-prime torsion degree"); torsion.push_back(f);
    }
    long period_lcm=1; std::vector<long> periods;
    for(long g=0;g<target.ell;++g) {
        Poly rem{one()}; long step;
        for(step=1;step<=period_cap;++step) {
            Poly shifted(source_M,E(0)); shifted.insert(shifted.end(),rem.begin(),rem.end());
            rem=right_divide(shifted,torsion[g],B).second;
            if(rem==Poly{one()}) break;
            if(step%128==0) std::cerr<<label<<": Frobenius group "<<g<<" step "<<step<<'\n';
        }
        require(step<=period_cap,"RESOURCE_LIMIT: Frobenius period search limit");
        periods.push_back(step); long gcd=std::gcd(period_lcm,step);
        require(period_lcm/gcd<=max_M/step,"RESOURCE_LIMIT: joint splitting degree exceeds max_M");
        period_lcm=period_lcm/gcd*step;
    }
    long joint=source_M*period_lcm;
    long M=joint*((std::max(target.M,target.n())+joint-1)/joint);
    require(!strict_M||M==target.M,"SPLITTING_MISMATCH: requested M="+std::to_string(target.M)+" is not a multiple of the joint splitting degree "+std::to_string(joint));
    require(M<=max_M,"RESOURCE_LIMIT: required extension degree "+std::to_string(M)+" exceeds max_M");
    for(auto& f:torsion) require(right_divide(psub(monomial(M),Poly{one()}),f,B).second.empty(),"splitting degree divisibility check");
    // Serialize before changing NTL's active extension-field modulus.
    std::vector<ZZ> source_mod;
    for(long j=0;j<=B.D;++j) source_mod.push_back(ZZ(NTL::rep(NTL::coeff(B.modulus,j))));
    auto archive=[&](const Poly& f) { std::vector<ZZ> z; for(auto& x:f) z.push_back(pack(x,B.p)); return z; };
    auto qw=pack(B.w,B.p); std::vector<std::vector<ZZ>> qsave,tsave;
    for(auto& f:Q) qsave.push_back(archive(f)); for(auto& f:torsion) tsave.push_back(archive(f));
    auto piSave=archive(pi),rhoSave=archive(rho),PIsave=archive(PI),Usave=archive(U),Vsave=archive(Vop),Wsave=archive(W);
    std::map<std::string,std::string> audit;
    audit["source_degree"]=std::to_string(source_M); audit["minimum_requested_M"]=std::to_string(target.M);
    audit["joint_splitting_degree"]=std::to_string(joint); audit["construction_seed"]=std::to_string(construction_seed);
    audit["theta_power"]=std::to_string(theta_power); audit["strict_M"]=strict_M?"1":"0";
    { std::ostringstream os; os<<pack_modulus(B.modulus,B.p); audit["source_modulus"]=os.str(); }
    { std::ostringstream os; for(auto v:periods) os<<v<<';'; audit["frobenius_periods"]=os.str(); }
    std::cerr<<label<<": target field M="<<M<<" (requested minimum "<<target.M<<")\n";
    Code C; target.M=M; C.p=target; C.family="hhy_automatic"; C.label=label;
    C.field=std::make_unique<Field>(target.q,M); Field& K=*C.field;
    NTL::zz_pEX source_polynomial;
    for(long j=0;j<static_cast<long>(source_mod.size());++j) NTL::SetCoeff(source_polynomial,j,E(NTL::conv<long>(source_mod[j])));
    NTL::vec_zz_pE roots; NTL::FindRoots(roots,source_polynomial);
    require(roots.length()==static_cast<long>(source_mod.size())-1,"source embedding root count");
    E root=roots[0]; for(long j=1;j<roots.length();++j) if(pack(roots[j],K.p)<pack(root,K.p)) root=roots[j];
    auto lift_element=[&](ZZ z) { E sum,term(1); while(z>0) { sum+=E(NTL::conv<long>(z%K.p))*term; z/=K.p; term*=root; } return sum; };
    auto lift=[&](const std::vector<ZZ>& f) { Poly v; for(auto& z:f) v.push_back(lift_element(z)); return v; };
    K.w=lift_element(qw); K.coordinates();
    Q.clear(); for(auto& f:qsave) Q.push_back(lift(f)); torsion.clear(); for(auto& f:tsave) torsion.push_back(lift(f));
    Poly piK=lift(piSave),rhoK=lift(rhoSave),PIK=lift(PIsave),UK=lift(Usave),VK=lift(Vsave),WK=lift(Wsave);
    for(long g=0;g<target.ell;++g) {
        std::cerr<<label<<": torsion kernel "<<g<<'\n';
        V beta=root_basis(torsion[g],K); require(sz(beta)==target.R,"torsion root-space dimension");
        long qi=g==0?target.s:0;
        E residue=g==0?one():oevaluate(piK,K.fq[g-1])/epow(oevaluate(rhoK,K.fq[g-1]),nu);
        require(!zero(residue),"zero auxiliary-function residue"); V gamma;
        for(auto& x:beta) {
            E y=evaluate(Q[qi],x,K)/residue;
            require(evaluate(PIK,y,K)==x,"inverse-isogeny point check");
            if(g==0) require(zero(evaluate(UK,y,K))&&zero(evaluate(VK,y,K)),"source infinity torsion");
            else {
                E z=K.fq[g-1]; require(zero(evaluate(psub(VK,pscale(z,UK)),y,K))&&zero(evaluate(psub(WK,pscale(z,VK)),y,K)),"source finite torsion");
            } gamma.push_back(y);
        }
        C.beta.insert(C.beta.end(),beta.begin(),beta.end()); C.gamma.push_back(gamma);
    }
    C.polys=message_cut(Q,target,K); C.audit=std::move(audit);
    C.audit["construction"]="HHY_Example_4.29_precomposition_degree_two_ring";
    C.audit["rng_note"]="C++_seeded_construction_not_bit_identical_to_Sage";
    C.finish(); return C;
}

// Standard polynomial control: W_i=z^i F_(q^R). Message monomials have
// q-exponents hR+j (h<s,j<r) and sR+j (j<b). This is not HHY arithmetic.
inline Code build_polynomial(Parameters p,const std::string& label) {
    p.check(); require(p.M%p.R==0&&p.ell<=p.M/p.R,"polynomial control needs M divisible by R and ell<=M/R");
    Code C; C.p=p; C.family="polynomial_control"; C.label=label; C.field=std::make_unique<Field>(p.q,p.M); auto& K=*C.field;
    NTL::zz_pX small; NTL::BuildIrred(small,K.a*p.R); NTL::zz_pEX lifted;
    for(long j=0;j<=NTL::deg(small);++j) NTL::SetCoeff(lifted,j,E(NTL::coeff(small,j)));
    NTL::vec_zz_pE roots; NTL::FindRoots(roots,lifted);
    require(roots.length()==NTL::deg(small),"polynomial subfield roots"); E alpha=roots[0];
    for(long j=1;j<roots.length();++j) if(pack(roots[j],K.p)<pack(alpha,K.p)) alpha=roots[j];
    V sub(p.R,one()); for(long j=1;j<p.R;++j) sub[j]=sub[j-1]*alpha;
    require(K.rankq(sub)==p.R,"polynomial subfield basis rank");
    for(long g=0;g<p.ell;++g) {
        V block; for(auto& a:sub) block.push_back(K.basis[g]*a);
        C.beta.insert(C.beta.end(),block.begin(),block.end()); C.gamma.push_back(sub);
    }
    for(long h=0;h<=p.s;++h) for(long j=0;j<(h==p.s?p.b:p.r);++j) C.polys.push_back(monomial(h*p.R+j));
    C.audit["construction"]="subfield_polynomial_rank_locality_control";
    C.audit["historical_note"]="same_parameter_control_reconstruction_not_claimed_identical_historical_points";
    C.finish(); return C;
}
} // namespace hhy
