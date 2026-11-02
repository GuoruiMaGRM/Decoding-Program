#pragma once
#include "decoding.hpp"
#include <NTL/RR.h>

namespace hhy {
inline uint64_t stable_hash(const std::string& s) {
    uint64_t h=14695981039346656037ULL; for(unsigned char c:s) { h^=c; h*=1099511628211ULL; } return h;
}
inline std::vector<long> sample_indices(long n,long t,std::mt19937_64& rng) {
    require(0<=t&&t<=n,"sample size"); std::vector<long> ix(n); std::iota(ix.begin(),ix.end(),0);
    std::shuffle(ix.begin(),ix.end(),rng); ix.resize(t); return ix;
}
inline V independent(const Field& K,long t,std::mt19937_64& rng,V initial={}) {
    require(t<=K.M,"error rank exceeds extension degree");
    while(sz(initial)<t) { E x=K.random(rng); V v=initial; v.push_back(x); if(K.rankq(v)==sz(v)) initial.push_back(x); } return initial;
}
inline V error_vector(const Code& C,long t,const std::string& kind,std::mt19937_64& rng) {
    const auto& K=*C.field; const auto& p=C.p; require(0<=t&&t<=std::min(p.n(),p.M),"error rank range");
    if(t==0) return V(p.n());
    if(kind=="column_sparse") {
        V basis=independent(K,t,rng),e(p.n()); auto ix=sample_indices(p.n(),t,rng);
        for(long j=0;j<t;++j) e[ix[j]]=basis[j]; return e;
    }
    if(kind=="crisscross") {
        long nr=t/2,nc=t-nr; V row;
        for(long j:sample_indices(K.M,nr,rng)) row.push_back(K.basis[j]);
        V basis=independent(K,t,rng,row); auto ix=sample_indices(p.n(),nc,rng);
        for(long attempt=0;attempt<100000;++attempt) {
            V e(p.n()); for(long j=0;j<p.n();++j) for(long h=0;h<nr;++h) e[j]+=basis[h]*K.randomq(rng);
            for(long h=0;h<nc;++h) e[ix[h]]+=basis[nr+h];
            if(K.rankq(e)==t) return e;
        } throw std::runtime_error("crisscross rank rejection limit");
    }
    V basis;
    if(kind=="row_sparse") for(long j:sample_indices(K.M,t,rng)) basis.push_back(K.basis[j]);
    else basis=independent(K,t,rng);
    std::vector<long> positions(p.n()); std::iota(positions.begin(),positions.end(),0);
    if(kind=="first_block"||kind=="last_block") {
        require(t<=p.R,"localized error rank exceeds R"); positions.resize(p.R);
        long start=kind=="first_block"?0:(p.ell-1)*p.R; for(long j=0;j<p.R;++j) positions[j]=start+j;
    }
    for(long attempt=0;attempt<100000;++attempt) {
        V e(p.n());
        for(long j:positions) for(long h=0;h<t;++h) e[j]+=basis[h]*K.randomq(rng);
        if(K.rankq(e)!=t) continue;
        if(kind=="support_captured"&&K.rankq(prefix(e,p.R))!=t) continue;
        return e;
    } throw std::runtime_error("fixed-rank error rejection limit");
}
struct Model { std::string kind; long t; bool erased=false; std::string name() const { return erased?kind:kind+"_t"+std::to_string(t); } };
inline std::vector<Model> models(const Parameters& p,const std::string& ranks,bool captured,bool erasures) {
    std::set<long> tt;
    if(ranks=="all") for(long t=0;t<=p.tau();++t) tt.insert(t);
    else if(ranks=="paper") tt={0,1,2,3,p.u(),p.u()+1,p.tau()/2,p.R,p.tau()};
    else {
        std::istringstream in(ranks); std::string s;
        while(std::getline(in,s,';')) if(!s.empty()) tt.insert(std::stol(s));
    }
    std::vector<Model> out;
    for(auto kind:{"uniform_rank","row_sparse","column_sparse","first_block","last_block","crisscross"})
        for(long t:tt) if(t>=0&&t<=p.tau()&&(t>0||std::string(kind)=="uniform_rank")&&
             (t<=p.R||(std::string(kind)!="first_block"&&std::string(kind)!="last_block"))) out.push_back({kind,t,false});
    if(captured) for(long t:tt) if(t>0&&t<=std::min({p.tau(),p.R,p.u()})) out.push_back({"support_captured",t,false});
    if(erasures) { out.push_back({"one_block_erasure",-1,true}); out.push_back({"all_blocks_erasure",-1,true}); }
    return out;
}
inline ZZ phi(long q,long a,long m) {
    if(a<0||m<0||a>m) return ZZ(0); ZZ z(1),qm=integer_power(q,m);
    for(long j=0;j<a;++j) z*=qm-integer_power(q,j); return z;
}
inline ZZ gaussian(long q,long m,long v) {
    if(v<0||v>m||m<0) return ZZ(0); return phi(q,v,m)/phi(q,v,v);
}
inline double exact_ratio(const ZZ& a,const ZZ& b) {
    require(b>0,"zero probability denominator"); return NTL::conv<double>(NTL::conv<NTL::RR>(a)/NTL::conv<NTL::RR>(b));
}
struct Theory { ZZ denominator,local_good,good_captured,captured; };
inline Theory theory(const Parameters& p,long t) {
    std::vector<ZZ> D(t+1); D[0]=1;
    for(long h=0;h<p.s;++h) {
        std::vector<ZZ> next(t+1);
        for(long v=0;v<=t;++v) if(D[v]!=0) for(long j=0;j<=std::min({p.u(),t,p.R});++j)
            for(long z=0;z<=std::min(v,j);++z) if(j-z<=t-v) {
                long w=v+j-z;
                ZZ count=gaussian(p.q,v,z)*gaussian(p.q,t-v,j-z)*integer_power(p.q,(v-z)*(j-z))*phi(p.q,j,p.R);
                next[w]+=D[v]*count;
            }
        D.swap(next);
    }
    long L=p.n()-p.s*p.R; Theory out; out.denominator=phi(p.q,t,p.n());
    for(long v=0;v<=t;++v) out.local_good+=D[v]*integer_power(p.q,v*L)*phi(p.q,t-v,L);
    out.good_captured=D[t]*integer_power(p.q,t*L);
    out.captured=phi(p.q,t,p.s*p.R)*integer_power(p.q,t*L); return out;
}
inline void write_theory(const Code& C,const std::vector<Model>& specs,const std::filesystem::path& directory) {
    std::ofstream out(directory/"theory.csv"),counts(directory/"theory_exact_counts.tsv");
    out<<"case,q,M,n,R,r,b,ell,s,tau,u,t,P_first_local_good,P_first_good_captured,P_first_captured\n";
    counts<<"case\tt\tdenominator\tlocal_good\tgood_captured\tcaptured\n";
    std::set<long> seen;
    for(auto& m:specs) if(!m.erased&&seen.insert(m.t).second) {
        auto a=theory(C.p,m.t); auto p=C.p;
        out<<std::setprecision(17)<<C.label<<','<<p.q<<','<<p.M<<','<<p.n()<<','<<p.R<<','<<p.r<<','<<p.b<<','<<p.ell<<','<<p.s<<','<<p.tau()<<','<<p.u()<<','<<m.t<<','
           <<exact_ratio(a.local_good,a.denominator)<<','<<exact_ratio(a.good_captured,a.denominator)<<','<<exact_ratio(a.captured,a.denominator)<<'\n';
        counts<<C.label<<'\t'<<m.t<<'\t'<<a.denominator<<'\t'<<a.local_good<<'\t'<<a.good_captured<<'\t'<<a.captured<<'\n';
    }
}
inline void progress(const std::filesystem::path& dir,const std::string& stage) {
    auto tmp=dir/"progress.txt.tmp"; { std::ofstream o(tmp); o<<stage<<'\n'; }
    std::filesystem::rename(tmp,dir/"progress.txt");
}
struct Input {
    V y,truth,error;
    long first_span=0;
    bool first_good=false,first_captured=false,any_good=false;
    std::vector<long> block_ranks;
    ErasurePlan local_erasure,global_erasure;
    std::set<long> erased;
};
inline Input make_input(const Code& C,const Decoder& decoder,const Model& m,std::mt19937_64& rng) {
    auto& K=*C.field; auto& p=C.p; Input x; V msg;
    for(long j=0;j<p.k();++j) msg.push_back(K.random(rng)); x.truth=C.encode(msg);
    if(m.erased) {
        for(long g=0;g<(m.kind=="one_block_erasure"?1:p.ell);++g)
            for(long j:sample_indices(p.R,p.R-p.r,rng)) x.erased.insert(g*p.R+j);
        x.y=x.truth; for(long j:x.erased) x.y[j]=0;
        x.local_erasure=decoder.erasure_plan(x.erased,true); x.global_erasure=decoder.erasure_plan(x.erased,false);
        return x;
    }
    x.error=error_vector(C,m.t,m.kind,rng); require(K.rankq(x.error)==m.t,"generated rank differs from requested rank");
    x.y=plus(x.truth,x.error);
    for(long g=0;g<p.ell;++g) x.block_ranks.push_back(K.rankq(V(x.error.begin()+g*p.R,x.error.begin()+(g+1)*p.R)));
    if(p.s==0) { x.first_good=true; x.first_captured=m.t==0; x.first_span=0; x.any_good=true; }
    else for(size_t i=0;i<decoder.sets.size();++i) {
        auto& S=decoder.sets[i]; bool good=true; for(long g:S.groups) good&=x.block_ranks[g]<=p.u();
        x.any_good|=good;
        if(i==0) { x.first_good=good; x.first_span=K.rankq(take(x.error,S.chosen)); x.first_captured=x.first_span==m.t; }
    }
    return x;
}
inline void validate_answer(const Code& C,const Model& m,const Input& x,const std::string& method,const Decoded& a) {
    if(method=="hybrid_no_fallback"&&a.word.empty()) return;
    require(a.word==x.truth,"decoder output mismatch: "+method+" / "+m.name());
    if(!m.erased) require(C.certify(x.y,a.word),"candidate certificate failed: "+method);
    if(method=="hybrid"&&!m.erased&&m.t>0&&C.p.s>0) {
        require(!(x.any_good&&a.stats.fallback),"sufficient local condition contradicted by fallback");
        require(!(x.first_good&&x.first_captured&&a.stats.branch!="support"),"full error-space capture did not return by interpolation");
    }
}

struct BenchOptions {
    uint64_t seed=1001;
    long inputs=16,rounds=3,batch=4,warmup=1,dense_cases=1;
    std::string ranks="all";
    bool ablations=true,captured=true,erasures=true,dump_inputs=true,gabidulin_controls=true;
};
inline std::string raw_header() {
    return "case,family,q,M,R,r,b,ell,s,n,k,KG,tau,u,model,t,seed,input,round,method,elapsed_ns,correct,branch,fallback,local_calls,attempts,first_local_good,first_captured,first_good_captured,first_span,any_local_good,conditional_model,degenerate_model";
}
inline void benchmark(const Code& C,const std::filesystem::path& dir,const BenchOptions& opt,DecoderOptions dop) {
    require(opt.inputs>0&&opt.rounds>0&&opt.batch>0&&opt.dense_cases>=0,"invalid benchmark counts");
    std::filesystem::create_directories(dir); progress(dir,"decoder_setup");
    auto setup_start=std::chrono::steady_clock::now(); Decoder D(C,dop);
    auto setup_ns=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-setup_start).count();
    auto specs=models(C.p,opt.ranks,opt.captured,opt.erasures);
    write_theory(C,specs,dir);
    { std::ofstream out(dir/"setup.txt"); out<<"NTL="<<NTL_VERSION<<"\ncompiler="<<__VERSION__<<"\nseed="<<opt.seed<<"\nsetup_ns="<<setup_ns<<"\ninputs_per_model="<<opt.inputs<<"\nrounds="<<opt.rounds<<"\nmax_sets="<<dop.max_sets<<"\ninfo_attempts="<<dop.info_attempts<<"\nmax_trial_radius="<<dop.max_trial_radius<<"\n"; }
    std::ofstream combined(dir/"raw.csv"); combined<<raw_header()<<'\n'; combined.flush();
    std::ofstream controls(dir/"gabidulin_controls.csv");
    controls<<"case,code_family,q,M,n,dimension,distance,common_test_radius,model,t,seed,input,round,method,elapsed_ns,correct\n";
    std::shared_ptr<ComponentPlan> rate_control;
    if(opt.gabidulin_controls) rate_control=std::make_shared<ComponentPlan>(C.beta,C.p.k(),C.p.tau(),*C.field);
    std::ofstream inputs;
    if(opt.dump_inputs) { inputs.open(dir/"inputs.txt"); inputs<<"# Packed polynomial-basis field elements. Private diagnostic sent words are not decoder inputs.\n"; }
    for(auto& model:specs) {
        progress(dir,"generate/"+model.name()); std::cerr<<C.label<<" seed="<<opt.seed<<" model="<<model.name()<<'\n';
        std::mt19937_64 rng(opt.seed^stable_hash(C.label+"/"+model.name())); std::vector<Input> samples;
        for(long i=0;i<opt.inputs;++i) {
            progress(dir,"generate/"+model.name()+"/"+std::to_string(i)); samples.push_back(make_input(C,D,model,rng));
            if(opt.dump_inputs) {
                inputs<<model.name()<<' '<<i<<"\ntruth "; write_vec(inputs,samples.back().truth,C.field->p);
                inputs<<"received "; write_vec(inputs,samples.back().y,C.field->p);
                inputs<<"error "; write_vec(inputs,samples.back().error,C.field->p);
                inputs<<"erased "<<samples.back().erased.size(); for(long j:samples.back().erased) inputs<<' '<<j; inputs<<'\n';
            }
        }
        std::vector<std::string> active=model.erased?std::vector<std::string>{"local_erasure","global_erasure"}:Decoder::methods(opt.ablations);
        auto call=[&](const std::string& method,const Input& x) {
            if(model.erased) return D.erase(x.y,method=="local_erasure"?x.local_erasure:x.global_erasure);
            return D.call(method,x.y);
        };
        // Complete validation precedes timing. A mismatch aborts the seed worker.
        for(long i=0;i<opt.inputs;++i) {
            for(auto& method:active) {
                progress(dir,"validate/"+model.name()+"/"+std::to_string(i)+"/"+method);
                validate_answer(C,model,samples[i],method,call(method,samples[i]));
            }
            if(!model.erased&&i<opt.dense_cases) {
                progress(dir,"dense_reference/"+model.name()+"/"+std::to_string(i));
                auto a=dense_reference(*D.full,samples[i].y,*C.field);
                require(a&&a->word==samples[i].truth,"independent dense reference disagrees");
                if(!D.sets.empty()) {
                    auto& S=D.sets[0]; V e=take(samples[i].error,S.chosen),basis=C.field->independent_basis(e);
                    Poly A=annihilator(basis,*C.field);
                    require(C.field->rankq(evaluations(A,samples[i].error,*C.field))==model.t-sz(basis),"projected error-rank identity failed");
                }
            }
        }
        for(long i=0;i<std::min(opt.warmup,opt.inputs);++i) for(auto& method:active) {
            progress(dir,"warmup/"+model.name()); validate_answer(C,model,samples[i],method,call(method,samples[i]));
        }
        auto path=dir/(model.name()+".csv.partial"); std::ofstream out(path); out<<raw_header()<<'\n';
        std::mt19937_64 order_rng(opt.seed^stable_hash("timing/"+C.label+"/"+model.name()));
        for(long round=0;round<opt.rounds;++round) {
            std::vector<long> order(opt.inputs); std::iota(order.begin(),order.end(),0); std::shuffle(order.begin(),order.end(),order_rng);
            for(long start=0;start<opt.inputs;start+=opt.batch) {
                std::shuffle(active.begin(),active.end(),order_rng);
                for(auto& method:active) {
                    progress(dir,"timing/"+model.name()+"/"+method+"/"+std::to_string(round)+"/"+std::to_string(start));
                    for(long pos=start;pos<std::min(start+opt.batch,opt.inputs);++pos) {
                        long i=order[pos]; auto before=std::chrono::steady_clock::now();
                        Decoded a=call(method,samples[i]);
                        auto elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-before).count();
                        validate_answer(C,model,samples[i],method,a); // truth comparison is outside timed region.
                        auto p=C.p; auto& x=samples[i]; bool correct=a.word==x.truth;
                        std::ostringstream row;
                        row<<C.label<<','<<C.family<<','<<p.q<<','<<p.M<<','<<p.R<<','<<p.r<<','<<p.b<<','<<p.ell<<','<<p.s<<','<<p.n()<<','<<p.k()<<','<<p.kg()<<','<<p.tau()<<','<<p.u()<<','
                           <<model.name()<<','<<model.t<<','<<opt.seed<<','<<i<<','<<round<<','<<method<<','<<elapsed<<','<<correct<<','<<a.stats.branch<<','<<a.stats.fallback<<','<<a.stats.local_calls<<','<<a.stats.attempts<<','
                           <<x.first_good<<','<<x.first_captured<<','<<(x.first_good&&x.first_captured)<<','<<x.first_span<<','<<x.any_good<<','<<(model.kind=="support_captured")<<','<<(model.kind=="crisscross"&&model.t==1)<<'\n';
                        out<<row.str();
                    } out.flush();
                }
            }
        }
        out.close(); auto final=dir/(model.name()+".csv"); std::filesystem::rename(path,final);
        { std::ifstream done(final); std::string line; std::getline(done,line); while(std::getline(done,line)) combined<<line<<'\n'; combined.flush(); }
        // Separate code comparison, never pooled with same-code decoder ratios.
        // Translate each identical error onto a fresh random Gabidulin codeword.
        if(rate_control&&!model.erased) {
            for(auto family:{"dimension_matched_gabidulin","ambient_gabidulin"}) {
                const auto& P=std::string(family)=="ambient_gabidulin"?*D.full:*rate_control;
                std::vector<V> words,received;
                std::mt19937_64 grng(opt.seed^stable_hash(C.label+"/"+model.name()+"/"+family));
                for(auto& x:samples) {
                    Poly f; for(long j=0;j<P.k;++j) f.push_back(C.field->random(grng));
                    words.push_back(evaluations(f,C.beta,*C.field)); received.push_back(plus(words.back(),x.error));
                }
                for(long round=0;round<opt.rounds;++round) for(long i=0;i<opt.inputs;++i) {
                    std::vector<std::string> gm{"incremental","gao_cpp"}; std::shuffle(gm.begin(),gm.end(),order_rng);
                    for(auto& method:gm) {
                        progress(dir,"separate_gabidulin_control/"+model.name()+"/"+family+"/"+method);
                        auto before=std::chrono::steady_clock::now();
                        auto a=method=="incremental"?incremental(P,received[i],*C.field):gao_cpp(P,received[i],*C.field);
                        auto ns=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-before).count();
                        require(a&&a->word==words[i],"separate Gabidulin control failed");
                        controls<<C.label<<','<<family<<','<<C.p.q<<','<<C.p.M<<','<<C.p.n()<<','<<P.k<<','<<C.p.n()-P.k+1<<','<<C.p.tau()<<','<<model.name()<<','<<model.t<<','<<opt.seed<<','<<i<<','<<round<<','<<method<<','<<ns<<",1\n";
                    }
                }
            } controls.flush();
        }
    }
    combined.close(); progress(dir,"complete"); { std::ofstream o(dir/"SUCCESS"); o<<"all requested models validated\n"; }
}
} // namespace hhy
