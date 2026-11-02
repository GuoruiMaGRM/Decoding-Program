#include "report.hpp"

using namespace hhy;
struct Arguments {
    std::map<std::string,std::string> kv;
    explicit Arguments(int argc,char** argv) {
        for(int i=1;i<argc;++i) {
            std::string key=argv[i]; require(key.rfind("--",0)==0,"expected --option");
            if(key=="--help") { kv[key]="1"; continue; }
            require(i+1<argc,"missing value for "+key); kv[key]=argv[++i];
        }
    }
    std::string get(const std::string& k,const std::string& fallback="") const { auto i=kv.find("--"+k); return i==kv.end()?fallback:i->second; }
    long number(const std::string& k,long fallback) const { return std::stol(get(k,std::to_string(fallback))); }
    bool flag(const std::string& k,bool fallback) const { auto s=get(k,fallback?"1":"0"); require(s=="0"||s=="1","boolean option "+k+" must be 0 or 1"); return s=="1"; }
};

static void verify(const Code& C) {
    const auto& K=*C.field; std::mt19937_64 rng(312479);
    for(long i=0;i<12;++i) {
        E x=K.random(rng),y=K.random(rng);
        require(K.from_coords(K.coords(x))==x,"coordinate roundtrip failed");
        if(!zero(x)) require(x*NTL::inv(x)==one(),"finite-field inverse failed");
        require(K.frob(x+y)==K.frob(x)+K.frob(y),"Frobenius additivity");
        Poly f,g; for(long j=0;j<3;++j) { f.push_back(K.random(rng)); g.push_back(K.random(rng)); }
        trim(f); trim(g);
        require(evaluate(compose(f,g,K),x,K)==evaluate(f,evaluate(g,x,K),K),"composition evaluation");
        auto fq=left_quotient(compose(f,g,K),f,K); require(fq&&*fq==g,"left composition quotient");
        require(exact_right(compose(f,g,K),g,K)==f,"right composition quotient");
    }
    if(K.a>1) {
        require(K.rankq(V{one(),K.w})==1,"rank was measured over the prime field instead of F_q");
        require(NTL::deg(NTL::rep(K.w))>0,"nonprime subfield generator degenerated");
    }
    Decoder D(C);
    for(long t:std::set<long>{0L,1L,C.p.u(),C.p.u()+1,C.p.tau()}) if(0<=t&&t<=C.p.tau()) {
        for(auto kind:{"uniform_rank","row_sparse","column_sparse","first_block","last_block","crisscross"}) {
            if((std::string(kind)=="first_block"||std::string(kind)=="last_block")&&t>C.p.R) continue;
            Model m{kind,t,false}; auto x=make_input(C,D,m,rng);
            auto ref=dense_reference(*D.full,x.y,K); require(ref&&ref->word==x.truth,"dense correctness reference failed");
            for(auto& method:Decoder::methods(true)) validate_answer(C,m,x,method,D.call(method,x.y));
        }
    }
    // Outside the guaranteed radius an accepted word must still be certified.
    // It need not be the original word, and failure is not forced.
    if(C.p.tau()+1<=C.p.n()) {
        V truth=C.encode(V(C.p.k(),one())); V e=error_vector(C,C.p.tau()+1,"uniform_rank",rng),y=plus(truth,e);
        for(auto& method:Decoder::methods(true)) { auto a=D.call(method,y); if(!a.word.empty()) require(C.certify(y,a.word),"uncertified output beyond radius"); }
    }
    std::cout<<"Verification passed for "<<C.label<<'\n';
}

int main(int argc,char** argv) {
    try {
        Arguments a(argc,argv);
        if(a.flag("help",false)) {
            std::cout<<"C++17/NTL HHY experiment suite\n"
                     <<"--task prepare --family hhy|polynomial --q Q --M MIN_M --R R --r r --b b --ell E --s S --source-M 2|4 --label ID --output FILE\n"
                     <<"--task verify --code FILE\n"
                     <<"--task bench --code FILE --output DIR --seed N --inputs 16 --rounds 3 --ranks all|paper|0;1;2 --ablations 1\n"
                     <<"--task aggregate --root DIR --expected-seeds 20 --bootstrap 2000\n"
                     <<"--task order --plan FILE --seed 8422\n"; return 0;
        }
        NTL::SetSeed(ZZ(a.number("construction-seed",20261027)));
        auto task=a.get("task");
        if(task=="prepare") {
            Parameters p; p.q=a.number("q",2); p.M=a.number("M",16); p.R=a.number("R",4); p.r=a.number("r",2);
            p.b=a.number("b",1); p.ell=a.number("ell",3); p.s=a.number("s",1);
            std::string label=a.get("label","case"),family=a.get("family","hhy");
            require(label.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")==std::string::npos,"label must be filename-safe");
            Code C;
            if(family=="hhy") {
                long limit=a.number("construction-attempts",16); require(limit>0,"construction attempt limit must be positive");
                bool done=false;
                for(long attempt=0;attempt<limit;++attempt) {
                    try {
                        C=build_hhy(p,a.number("source-M",p.q==7?2:4),a.number("max-M",512),a.number("period-cap",4096),
                                    a.number("construction-seed",20261027)+1009*attempt,label,a.flag("strict-M",true),1+attempt);
                        C.audit["source_search_attempts"]=std::to_string(attempt+1); done=true; break;
                    } catch(const std::runtime_error& e) {
                        std::string message=e.what();
                        bool selector=message.rfind("SPLITTING_MISMATCH:",0)==0||message.rfind("SOURCE_SELECTOR:",0)==0;
                        if(!selector) throw;
                        std::cerr<<label<<": source selection "<<attempt+1<<": "<<message<<'\n';
                    }
                }
                require(done,"SOURCE_SELECTION_LIMIT: no compatible source found; no code or timing was substituted");
            }
            else if(family=="polynomial") C=build_polynomial(p,label);
            else throw std::runtime_error("unknown construction family");
            auto path=std::filesystem::path(a.get("output")); require(!path.empty(),"missing --output");
            if(!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
            auto tmp=path; tmp+=".partial"; save_code(C,tmp); std::filesystem::rename(tmp,path);
            std::cout<<"Prepared "<<C.family<<" "<<C.label<<" q="<<C.p.q<<" M="<<C.p.M<<" n="<<C.p.n()<<" k="<<C.p.k()<<" KG="<<C.p.kg()<<'\n';
        } else if(task=="verify") {
            Code C=load_code(a.get("code")); verify(C);
        } else if(task=="bench") {
            Code C=load_code(a.get("code")); BenchOptions o;
            o.seed=a.number("seed",1001); o.inputs=a.number("inputs",16); o.rounds=a.number("rounds",3); o.batch=a.number("batch",4);
            o.warmup=a.number("warmup",1); o.dense_cases=a.number("dense-cases",1); o.ranks=a.get("ranks","all");
            o.ablations=a.flag("ablations",true); o.captured=a.flag("captured",true); o.erasures=a.flag("erasures",true);
            o.dump_inputs=a.flag("dump-inputs",true); o.gabidulin_controls=a.flag("gabidulin-controls",true);
            DecoderOptions d; d.max_sets=a.number("max-sets",3); d.info_attempts=a.number("info-attempts",7); d.max_trial_radius=a.number("max-trial-radius",3);
            d.cyclic=a.get("candidate-policy","cyclic")=="cyclic";
            require(a.get("candidate-policy","cyclic")=="cyclic"||a.get("candidate-policy")=="lexicographic","unknown candidate policy");
            benchmark(C,a.get("output"),o,d);
        } else if(task=="aggregate") {
            aggregate(a.get("root"),a.number("expected-seeds",20),a.number("bootstrap",2000));
        } else if(task=="order") {
            std::ifstream in(a.get("plan")); require(in.good(),"cannot read plan");
            std::string line; std::vector<std::string> lines;
            while(std::getline(in,line)) if(!line.empty()&&line[0]!='#') lines.push_back(line);
            std::mt19937_64 rng(a.number("seed",8422)); std::shuffle(lines.begin(),lines.end(),rng);
            for(auto& s:lines) std::cout<<s<<'\n';
        } else throw std::runtime_error("unknown or missing --task; use --help");
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"ERROR: "<<e.what()<<'\n'; return 1;
    }
}
