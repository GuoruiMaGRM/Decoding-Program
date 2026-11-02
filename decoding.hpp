#pragma once
#include "construction.hpp"

namespace hhy {
struct ComponentPlan {
    V x;
    long k,radius;
    Mat clear_inverse;
    std::vector<Poly> newton;
    V inverse_discrepancies;
    Poly all_annihilator;
    ComponentPlan(const V& points,long dimension,long rho,const Field& K):x(points),k(dimension),radius(rho) {
        require(1<=k&&k<=sz(x)&&0<=2*radius&&2*radius<=sz(x)-k,"invalid component parameters");
        clear_inverse=inverse(moore(prefix(x,k),k,K));
        Poly a{one()};
        for(auto& v:x) {
            E delta=evaluate(a,v,K); require(!zero(delta),"dependent component points");
            newton.push_back(a); inverse_discrepancies.push_back(NTL::inv(delta));
            a=compose(Poly{-epow(delta,K.q-1),one()},a,K);
        } all_annihilator=a;
    }
};
struct ComponentAnswer { Poly f; V word; };
inline std::optional<ComponentAnswer> component_recover(const Poly& N,const Poly& Vop,
                                                       const ComponentPlan& P,const V& y,const Field& K) {
    if(Vop.empty()||degree(Vop)>P.radius||degree(N)>=P.k+P.radius) return {};
    auto f=left_quotient(N,Vop,K);
    if(!f||degree(*f)>=P.k) return {};
    V word=evaluations(*f,P.x,K);
    if(K.rankq(minus(y,word))>P.radius) return {};
    return ComponentAnswer{*f,std::move(word)};
}

// Incremental two-row linearized interpolation (common reconstruction core).
// No true error rank, support, sent word, or error positions are passed here.
inline std::optional<ComponentAnswer> incremental(const ComponentPlan& P,const V& y,const Field& K,bool screen=false) {
    require(sz(y)==sz(P.x),"component received length");
    if(screen) {
        Poly f=mv(P.clear_inverse,prefix(y,P.k)); V word=evaluations(f,P.x,K);
        if(word==y) return ComponentAnswer{f,word};
    }
    using Row=std::pair<Poly,Poly>;
    std::vector<Row> row{{Poly{one()},{}},{ {},Poly{one()} }};
    auto key=[&](const Row& r) {
        long dN=r.first.empty()?-100000000L:degree(r.first);
        long dV=r.second.empty()?-100000000L:degree(r.second)+P.k-1;
        return std::make_pair(std::max(dN,dV),dV>=dN?1:0);
    };
    for(long j=0;j<sz(P.x);++j) {
        E d[2]={evaluate(row[0].first,P.x[j],K)-evaluate(row[0].second,y[j],K),
                evaluate(row[1].first,P.x[j],K)-evaluate(row[1].second,y[j],K)};
        if(zero(d[0])&&zero(d[1])) continue;
        long pivot=zero(d[0])?1:(zero(d[1])?0:(key(row[0])<key(row[1])?0:1));
        Row old=row[pivot]; long other=1-pivot;
        if(!zero(d[other])) {
            E c=d[other]/d[pivot];
            row[other].first=psub(row[other].first,pscale(c,old.first));
            row[other].second=psub(row[other].second,pscale(c,old.second));
        }
        Poly factor{-epow(d[pivot],K.q-1),one()};
        row[pivot]={compose(factor,old.first,K),compose(factor,old.second,K)};
    }
    if(key(row[1])<key(row[0])) std::swap(row[0],row[1]);
    for(auto& r:row) { auto a=component_recover(r.first,r.second,P,y,K); if(a) return a; }
    return {};
}

// Independently written Gao-style decoder. Interpolate R, then perform the
// right Euclidean algorithm A=Q o R+rem, tracking LEFT multipliers of R.
// It is not Sage code and is not an optimized paper-author implementation.
inline std::optional<ComponentAnswer> gao_cpp(const ComponentPlan& P,const V& y,const Field& K) {
    require(sz(y)==sz(P.x),"Gao received length"); Poly R;
    for(long i=0;i<sz(y);++i) {
        E c=(y[i]-evaluate(R,P.x[i],K))*P.inverse_discrepancies[i];
        R=padd(R,pscale(c,P.newton[i]));
    }
    Poly r0=P.all_annihilator,r1=R,v0,v1{one()};
    while(!r1.empty()&&degree(r1)>=P.k+P.radius) {
        auto qr=right_divide(r0,r1,K);
        Poly v2=psub(v0,compose(qr.first,v1,K));
        r0=std::move(r1); r1=std::move(qr.second); v0=std::move(v1); v1=std::move(v2);
    }
    return component_recover(r1,v1,P,y,K);
}

// Independent dense homogeneous Welch--Berlekamp reference, excluded from
// reported timing. All nullspace vectors satisfying the degree bounds recover
// the same sent polynomial when the input is inside the decoding radius.
inline std::optional<ComponentAnswer> dense_reference(const ComponentPlan& P,const V& y,const Field& K) {
    Mat A;
    for(long i=0;i<sz(y);++i) {
        V row; E v=y[i];
        for(long j=0;j<=P.radius;++j) { row.push_back(-v); v=K.frob(v); }
        v=P.x[i]; for(long j=0;j<P.k+P.radius;++j) { row.push_back(v); v=K.frob(v); }
        A.push_back(row);
    }
    for(auto& z:kernel(A,P.k+2*P.radius+1)) {
        Poly Vop(z.begin(),z.begin()+P.radius+1),N(z.begin()+P.radius+1,z.end()); trim(Vop); trim(N);
        auto a=component_recover(N,Vop,P,y,K); if(a) return a;
    } return {};
}

struct Stats {
    std::string branch="failure";
    long local_calls=0,attempts=0;
    bool fallback=false;
};
struct Decoded { V word; Stats stats; };
struct DecoderOptions { long max_sets=3,info_attempts=7,max_trial_radius=3; bool cyclic=true; };
struct Shortening {
    std::vector<long> groups,chosen,keep;
    V eta;
    Mat transfer;
    std::map<long,Mat> support_inverse;
    std::map<long,std::shared_ptr<ComponentPlan>> projected;
    std::shared_ptr<ComponentPlan> unprojected;
};
struct ErasureMap { std::vector<long> input,output; Mat transfer; };
using ErasurePlan=std::vector<ErasureMap>;

class Decoder {
public:
    const Code& C;
    const Field& K;
    const Parameters& p;
    DecoderOptions config;
    std::shared_ptr<ComponentPlan> full,screened,small;
    std::vector<std::shared_ptr<ComponentPlan>> locals,trials;
    std::vector<std::pair<std::vector<long>,Mat>> information;
    std::vector<Shortening> sets;
    Mat global_prefix_inverse;
    explicit Decoder(const Code& c,DecoderOptions opt={}):C(c),K(*c.field),p(c.p),config(opt) {
        full=std::make_shared<ComponentPlan>(C.beta,p.kg(),p.tau(),K);
        screened=std::make_shared<ComponentPlan>(prefix(C.beta,p.kg()+2*p.tau()),p.kg(),p.tau(),K);
        long budget=std::min(p.tau(),std::max(1L,p.u()));
        if(budget<p.tau()) small=std::make_shared<ComponentPlan>(prefix(C.beta,p.kg()+2*budget),p.kg(),budget,K);
        for(long rho=1;rho<std::min(p.tau(),config.max_trial_radius+1);++rho)
            trials.push_back(std::make_shared<ComponentPlan>(prefix(C.beta,p.kg()+2*rho),p.kg(),rho,K));
        for(auto& x:C.gamma) locals.push_back(std::make_shared<ComponentPlan>(x,p.r,p.u(),K));
        global_prefix_inverse=inverse(moore(prefix(C.beta,p.kg()),p.kg(),K));
        for(long shift=0;shift<std::min(p.ell,config.info_attempts);++shift) {
            std::vector<long> order(p.n()); for(long j=0;j<p.n();++j) order[j]=(shift*p.R+j)%p.n();
            auto piv=rref(columns(C.G,order)).piv; require(static_cast<long>(piv.size())==p.k(),"information set rank");
            std::vector<long> ix; for(long j:piv) ix.push_back(order[j]);
            information.push_back({ix,inverse(transpose(columns(C.G,ix)))});
        }
        if(p.s==0) return;
        require(config.max_sets>0,"nonempty selected-block attempt list required");
        std::vector<std::vector<long>> selected;
        if(config.cyclic) {
            for(long shift=0;shift<std::min(p.ell,config.max_sets);++shift) {
                std::vector<long> g; for(long i=0;i<p.s;++i) g.push_back((shift+i)%p.ell); selected.push_back(g);
            }
        } else {
            std::vector<long> current;
            std::function<void(long)> enumerate=[&](long start) {
                if(static_cast<long>(selected.size())>=config.max_sets) return;
                if(static_cast<long>(current.size())==p.s) { selected.push_back(current); return; }
                for(long j=start;j<p.ell;++j) { current.push_back(j); enumerate(j+1); current.pop_back(); if(static_cast<long>(selected.size())>=config.max_sets) break; }
            }; enumerate(0);
        }
        for(auto& groups:selected) {
            Shortening S; S.groups=groups; Mat B;
            for(long g:groups) { for(long j=0;j<p.R;++j) S.chosen.push_back(g*p.R+j); B.insert(B.end(),C.local[g].begin(),C.local[g].end()); }
            std::set<long> gone(S.chosen.begin(),S.chosen.end());
            for(long j=0;j<p.n();++j) if(!gone.count(j)) S.keep.push_back(j);
            auto piv=rref(B).piv; require(static_cast<long>(piv.size())==p.s*p.r,"restriction map not onto");
            Mat inv=inverse(columns(B,piv)),section(p.k(),V(p.s*p.r));
            for(long i=0;i<static_cast<long>(piv.size());++i) section[piv[i]]=inv[i];
            S.transfer=mm(transpose(C.G),section);
            Poly A=annihilator(take(C.beta,S.chosen),K);
            S.eta=evaluations(A,take(C.beta,S.keep),K);
            require(K.rankq(S.eta)==sz(S.eta),"shortened evaluation points dependent");
            Mat kernel_words=columns(mm(kernel(B,p.k()),C.G),S.keep);
            Mat Gab=transpose(moore(S.eta,p.b,K)),joined=Gab;
            joined.insert(joined.end(),kernel_words.begin(),kernel_words.end());
            require(rows(kernel_words)==p.b&&rank(kernel_words)==p.b&&rank(joined)==p.b,"exact shortening identity failed");
            S.unprojected=std::make_shared<ComponentPlan>(prefix(S.eta,p.b+2*p.tau()),p.b,p.tau(),K);
            for(long a=0;a<=std::min(p.tau(),p.s*p.u());++a) {
                S.support_inverse[a]=inverse(moore(prefix(S.eta,p.b+a),p.b+a,K));
                S.projected[a]=std::make_shared<ComponentPlan>(prefix(S.eta,p.b+2*p.tau()-a),p.b+a,p.tau()-a,K);
            } sets.push_back(std::move(S));
        }
    }
    Decoded global(const V& y,bool optimized=true,bool clear=true,bool gao=false) const {
        if(clear&&C.member(y)) return {y,{"zero",0,0,false}};
        const ComponentPlan& P=*(optimized?screened:full);
        auto a=gao?gao_cpp(P,prefix(y,sz(P.x)),K):incremental(P,prefix(y,sz(P.x)),K);
        if(a) { V word=evaluations(a->f,C.beta,K); if(C.certify(y,word)) return {word,{gao?"gao":"global",0,0,false}}; }
        return {{},{"failure",0,0,false}};
    }
    Decoded fallback(const V& y,Stats stats) const {
        auto a=global(y,true,false); stats.fallback=true; stats.branch=a.word.empty()?"fallback_failure":"fallback";
        a.stats=stats; return a;
    }
    Decoded adaptive(const V& y) const {
        if(C.member(y)) return {y,{"zero",0,0,false}};
        if(small) {
            auto a=incremental(*small,prefix(y,sz(small->x)),K);
            if(a) { auto word=evaluations(a->f,C.beta,K); if(C.certify(y,word)) return {word,{"small_rank",0,1,false}}; }
        }
        if(!small) return global(y,true,false); return fallback(y,{"",0,1,false});
    }
    Decoded hybrid(const V& y,bool interpolate=true,bool project=true,long attempt_limit=-1,bool allow_fallback=true) const {
        if(C.member(y)) return {y,{"zero",0,0,false}};
        if(p.s==0) return global(y,true,false);
        Stats stats;
        long count=attempt_limit<0?static_cast<long>(sets.size()):std::min(attempt_limit,static_cast<long>(sets.size()));
        for(long si=0;si<count;++si) {
            ++stats.attempts; auto& S=sets[si]; V coeff,errors; bool good=true;
            for(long g:S.groups) {
                ++stats.local_calls; V part(y.begin()+g*p.R,y.begin()+(g+1)*p.R);
                auto a=incremental(*locals[g],part,K,true);
                if(!a) { good=false; break; }
                a->f.resize(p.r); coeff.insert(coeff.end(),a->f.begin(),a->f.end());
                auto e=minus(part,a->word); errors.insert(errors.end(),e.begin(),e.end());
            }
            if(!good) continue;
            V f0=mv(S.transfer,coeff),residual=minus(take(y,S.keep),take(f0,S.keep));
            auto combine=[&](const Poly& h) {
                V word=f0,eval=evaluations(h,S.eta,K); for(long j=0;j<sz(eval);++j) word[S.keep[j]]+=eval[j]; return word;
            };
            if(!project) {
                auto a=incremental(*S.unprojected,prefix(residual,sz(S.unprojected->x)),K);
                if(a) { V word=combine(a->f); if(C.certify(y,word)) { stats.branch="shortened_plain"; return {word,stats}; } }
                continue;
            }
            V basis=K.independent_basis(errors); long a=sz(basis);
            if(a>p.tau()||!S.support_inverse.count(a)) continue;
            Poly Lambda=annihilator(basis,K);
            if(interpolate) {
                V values=evaluations(Lambda,prefix(residual,p.b+a),K);
                Poly N=mv(S.support_inverse.at(a),values); trim(N);
                auto h=left_quotient(N,Lambda,K);
                if(h&&degree(*h)<p.b) {
                    V word=combine(*h); bool killed=true;
                    for(long j=0;j<p.n();++j) if(!zero(evaluate(Lambda,y[j]-word[j],K))) { killed=false; break; }
                    if(killed&&C.certify(y,word)) { stats.branch="support"; return {word,stats}; }
                }
            }
            auto& P=*S.projected.at(a);
            V transformed=evaluations(Lambda,prefix(residual,sz(P.x)),K);
            auto a2=incremental(P,transformed,K);
            if(a2) {
                auto h=left_quotient(a2->f,Lambda,K);
                if(h&&degree(*h)<p.b) { V word=combine(*h); if(C.certify(y,word)) { stats.branch="shortened"; return {word,stats}; } }
            }
        }
        if(allow_fallback) return fallback(y,stats);
        stats.branch="no_fallback_failure"; return {{},stats};
    }
    Decoded local_all(const V& y) const {
        V word; Stats stats;
        for(long g=0;g<p.ell;++g) {
            ++stats.local_calls; V part(y.begin()+g*p.R,y.begin()+(g+1)*p.R);
            auto a=incremental(*locals[g],part,K,true); if(!a) return fallback(y,stats);
            word.insert(word.end(),a->word.begin(),a->word.end());
        }
        if(C.certify(y,word)) { stats.branch="local_all"; return {word,stats}; } return fallback(y,stats);
    }
    Decoded information_sets(const V& y) const {
        Stats stats;
        for(auto& info:information) {
            ++stats.attempts; V word=C.encode(mv(info.second,take(y,info.first)));
            if(C.certify(y,word)) { stats.branch="information_set"; return {word,stats}; }
        } return fallback(y,stats);
    }
    Decoded prefix_interpolation(const V& y) const {
        auto f=mv(global_prefix_inverse,prefix(y,p.kg())); auto word=evaluations(f,C.beta,K);
        if(C.certify(y,word)) return {word,{"prefix",0,1,false}};
        return fallback(y,{"",0,1,false});
    }
    Decoded multi_radius(const V& y) const {
        if(C.member(y)) return {y,{"zero",0,0,false}};
        Stats stats;
        for(auto& P:trials) {
            ++stats.attempts; auto a=incremental(*P,prefix(y,sz(P->x)),K);
            if(a) { V word=evaluations(a->f,C.beta,K); if(C.certify(y,word)) { stats.branch="small_rank"; return {word,stats}; } }
        } return fallback(y,stats);
    }
    Decoded call(const std::string& method,const V& y) const {
        if(method=="hybrid") return hybrid(y);
        if(method=="global_full") return global(y,false,false);
        if(method=="global_screened") return global(y,true,true);
        if(method=="rank_adaptive") return adaptive(y);
        if(method=="local_all") return local_all(y);
        if(method=="information_sets") return information_sets(y);
        if(method=="prefix_interpolation") return prefix_interpolation(y);
        if(method=="multi_radius") return multi_radius(y);
        if(method=="gao_cpp") return global(y,false,false,true);
        if(method=="hybrid_no_interpolation") return hybrid(y,false);
        if(method=="hybrid_no_projection") return hybrid(y,false,false);
        if(method=="hybrid_one_set") return hybrid(y,true,true,1);
        if(method=="hybrid_no_fallback") return hybrid(y,true,true,-1,false);
        throw std::runtime_error("unknown decoder: "+method);
    }
    static std::vector<std::string> methods(bool ablations) {
        std::vector<std::string> m{"hybrid","global_full","global_screened","rank_adaptive","local_all","information_sets","prefix_interpolation","multi_radius","gao_cpp"};
        if(ablations) for(auto x:{"hybrid_no_interpolation","hybrid_no_projection","hybrid_one_set","hybrid_no_fallback"}) m.emplace_back(x);
        return m;
    }
    ErasurePlan erasure_plan(const std::set<long>& erased,bool local) const {
        ErasurePlan plans;
        for(long gg=0;gg<(local?p.ell:1);++gg) {
            std::vector<long> available,gone,relative;
            long start=local?gg*p.R:0,end=local?(gg+1)*p.R:p.n();
            for(long j=start;j<end;++j) if(erased.count(j)) gone.push_back(j); else { available.push_back(j); relative.push_back(j-start); }
            if(gone.empty()) continue;
            Mat G=local?transpose(moore(C.gamma[gg],p.r,K)):C.G;
            auto piv=rref(columns(G,relative)).piv; long dim=rows(G);
            require(static_cast<long>(piv.size())==dim,"unrecoverable erasure pattern");
            std::vector<long> inputs,local_input,local_output;
            for(long j:piv) { inputs.push_back(available[j]); local_input.push_back(available[j]-start); }
            for(long j:gone) local_output.push_back(j-start);
            Mat transfer=mm(transpose(columns(G,local_output)),inverse(transpose(columns(G,local_input))));
            plans.push_back({inputs,gone,transfer});
        } return plans;
    }
    Decoded erase(const V& y,const ErasurePlan& plans) const {
        V out=y; for(auto& P:plans) {
            auto v=mv(P.transfer,take(y,P.input)); for(size_t j=0;j<P.output.size();++j) out[P.output[j]]=v[j];
        } return {out,{"erasure",0,0,false}};
    }
};
} // namespace hhy
