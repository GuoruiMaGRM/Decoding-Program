#pragma once
#include "experiments.hpp"

namespace hhy {
inline std::vector<std::string> split(const std::string& s,char delimiter=',') {
    std::vector<std::string> a; std::string item; std::istringstream in(s);
    while(std::getline(in,item,delimiter)) { if(!item.empty()&&item.back()=='\r') item.pop_back(); a.push_back(item); } return a;
}
inline std::string escaped(const std::string& s) {
    std::string z; for(char c:s) { if(c=='&') z+="&amp;"; else if(c=='<') z+="&lt;"; else if(c=='>') z+="&gt;"; else if(c=='\"') z+="&quot;"; else z+=c; } return z;
}
struct Timings {
    std::vector<double> us;
    std::map<uint64_t,std::pair<long double,long>> process;
    std::set<std::pair<uint64_t,long>> distinct;
    long calls=0,correct=0,fall=0,support=0,shortened=0,zero_branch=0,local=0,first_good=0,first_captured=0,first_both=0;
    void add(uint64_t seed,long input,double ns,const std::vector<std::string>& r) {
        us.push_back(ns/1000.0); auto& p=process[seed]; p.first+=ns/1000.0; ++p.second;
        distinct.insert({seed,input}); ++calls; correct+=std::stol(r[21]); fall+=std::stol(r[23]);
        support+=r[22]=="support"; shortened+=r[22]=="shortened"; zero_branch+=r[22]=="zero";
        local+=std::stol(r[24]); first_good+=std::stol(r[26]); first_captured+=std::stol(r[27]); first_both+=std::stol(r[28]);
    }
    long double mean(uint64_t seed) const { auto p=process.at(seed); return p.first/p.second; }
    double mean() const { long double z=0; for(auto& p:process) z+=mean(p.first); return static_cast<double>(z/process.size()); }
    double fraction(long n) const { return calls?double(n)/calls:0; }
};
inline double quantile(std::vector<double> a,double q) {
    require(!a.empty(),"quantile of empty data"); std::sort(a.begin(),a.end()); double x=q*(a.size()-1); size_t i=static_cast<size_t>(x);
    return i+1==a.size()?a[i]:a[i]+(x-i)*(a[i+1]-a[i]);
}
struct Group {
    std::vector<std::string> meta;
    std::map<std::string,Timings> methods;
    bool all_workers_complete=true;
};
inline std::pair<double,double> bootstrap(const Group& g,const std::string& target,const std::vector<std::string>& baseline,
                                        const std::vector<uint64_t>& seeds,long resamples,uint64_t key) {
    if(seeds.size()<2||resamples<2) return {NAN,NAN};
    std::mt19937_64 rng(key); std::uniform_int_distribution<size_t> pick(0,seeds.size()-1); std::vector<double> values;
    for(long b=0;b<resamples;++b) {
        std::vector<uint64_t> chosen; for(size_t i=0;i<seeds.size();++i) chosen.push_back(seeds[pick(rng)]);
        long double denominator=0; for(auto s:chosen) denominator+=g.methods.at(target).mean(s);
        long double best=std::numeric_limits<long double>::infinity();
        for(auto& method:baseline) { long double num=0; for(auto s:chosen) num+=g.methods.at(method).mean(s); best=std::min(best,num); }
        values.push_back(static_cast<double>(best/denominator));
    } return {quantile(values,.025),quantile(values,.975)};
}
inline void aggregate(const std::filesystem::path& root,long expected_seeds,long resamples) {
    std::map<std::pair<std::string,std::string>,Group> groups;
    auto outdir=root/"report"; std::filesystem::create_directories(outdir);
    std::ofstream sources(outdir/"source_workers.tsv"); sources<<"file\tworker_complete\n";
    std::set<std::string> seen_sources;
    for(auto& entry:std::filesystem::recursive_directory_iterator(root)) if(entry.is_regular_file()&&entry.path().extension()==".csv"&&entry.path().filename()!="raw.csv"&&entry.path().parent_path().filename().string().rfind("seed_",0)==0) {
        auto path=entry.path(); std::ifstream in(path); std::string line; std::getline(in,line);
        if(line!=raw_header()) continue; // theory and separate-code controls have other schemas.
        sources<<path.string()<<'\t'<<std::filesystem::exists(path.parent_path()/"SUCCESS")<<'\n';
        while(std::getline(in,line)) {
            auto r=split(line); require(r.size()==33,"malformed raw row: "+path.string());
            std::string unique=r[0]+"/"+r[14]+"/"+r[16]+"/"+r[17]+"/"+r[18]+"/"+r[19];
            require(seen_sources.insert(unique).second,"duplicate timing identity; separate run directories must not be merged");
            auto& g=groups[{r[0],r[14]}]; if(g.meta.empty()) g.meta=r;
            g.all_workers_complete&=std::filesystem::exists(path.parent_path()/"SUCCESS");
            g.methods[r[19]].add(std::stoull(r[16]),std::stol(r[17]),std::stod(r[20]),r);
        }
    }
    require(!groups.empty(),"no completed model rows found");
    std::ofstream summary(outdir/"summary.csv"),methods(outdir/"methods.csv"),ratios(outdir/"ratios.csv"),html(outdir/"report.html");
    summary<<"case,family,q,M,R,r,b,ell,s,n,k,KG,tau,u,model,t,target,best_baseline,target_us,best_baseline_us,min_ratio,ci95_low,ci95_high,process_seeds,distinct_inputs,complete_seed_count,all_contributing_workers_complete,observed_color,claim_scope,conditional_model,degenerate_model\n";
    methods<<"case,model,method,process_mean_us,median_call_us,p95_call_us,process_seeds,distinct_inputs,timed_calls,success_fraction,fallback_fraction,support_fraction,shortened_fraction,zero_fraction,mean_local_calls,first_local_good_fraction,first_captured_fraction,first_good_captured_fraction,is_ablation\n";
    ratios<<"case,model,target,comparator,comparator_us,target_us,ratio,ci95_low,ci95_high,color,paired_process_seeds,is_ablation\n";
    html<<"<!doctype html><meta charset='utf-8'><title>HHY C++ NTL results</title><style>body{font-family:sans-serif;margin:24px}table{border-collapse:collapse;margin:16px 0}th,td{border:1px solid #999;padding:5px;text-align:right}th:first-child,td:first-child{text-align:left}.green{background:#cdebd1}.gray{background:#dedede}.missing{background:white}p{max-width:1100px}</style><h1>HHY C++/NTL results</h1><p>Green: comparator mean / hybrid mean &gt; 1. Gray: ratio &le; 1. Colors show the observed inequality, not statistical significance. The minimum uses the eight same-code baselines only. Ablations and different-code Gabidulin controls are excluded. Confidence intervals resample paired process seeds and reselect the fastest comparator. Repeated timings are not independent inputs. Known erasures and conditional support-captured errors remain separate.</p>";
    const std::vector<std::string> core{"global_full","global_screened","rank_adaptive","local_all","information_sets","prefix_interpolation","multi_radius","gao_cpp"};
    std::map<std::string,std::vector<std::pair<std::string,std::vector<double>>>> heat;
    for(auto& kv:groups) {
        auto& g=kv.second; auto& m=g.meta; bool erasure=m[15]=="-1"; std::string target=erasure?"local_erasure":"hybrid";
        auto baseline=erasure?std::vector<std::string>{"global_erasure"}:core;
        require(g.methods.count(target),"missing target method");
        std::vector<uint64_t> paired;
        for(auto& s:g.methods.at(target).process) {
            bool all=true; for(auto& b:baseline) all&=g.methods.count(b)&&g.methods.at(b).process.count(s.first);
            if(all) paired.push_back(s.first);
        }
        require(!paired.empty(),"no paired process seed for "+kv.first.first+"/"+kv.first.second);
        auto mean=[&](const std::string& method) { long double z=0; for(auto s:paired) z+=g.methods.at(method).mean(s); return static_cast<double>(z/paired.size()); };
        double tm=mean(target),best=std::numeric_limits<double>::infinity(); std::string best_name;
        require(tm>0,"nonpositive target timing");
        for(auto& b:baseline) if(mean(b)<best) { best=mean(b); best_name=b; }
        auto ci=bootstrap(g,target,baseline,paired,resamples,stable_hash(kv.first.first+kv.first.second));
        bool complete=static_cast<long>(paired.size())==expected_seeds;
        bool sufficient=complete&&g.all_workers_complete&&paired.size()>=20;
        std::string claim=sufficient?"paired_seed_exploratory":"incomplete_or_pilot";
        summary<<std::setprecision(17);
        for(long i=0;i<16;++i) summary<<m[i]<<',';
        summary<<target<<','<<best_name<<','<<tm<<','<<best<<','<<best/tm<<',';
        if(std::isfinite(ci.first)) summary<<ci.first; summary<<','; if(std::isfinite(ci.second)) summary<<ci.second;
        summary<<','<<paired.size()<<','<<g.methods.at(target).distinct.size()<<','<<complete<<','<<g.all_workers_complete<<','<<(best>tm?"green":"gray")<<','<<claim<<','<<m[31]<<','<<m[32]<<'\n';
        html<<"<h2>"<<escaped(m[0]+" / "+m[14])<<"</h2><p>q="<<m[2]<<", M="<<m[3]<<", R="<<m[4]<<", r="<<m[5]<<", b="<<m[6]<<", ell="<<m[7]<<", s="<<m[8]<<", n="<<m[9]<<", k="<<m[10]<<", KG="<<m[11]<<", tau="<<m[12]<<", u="<<m[13]<<"; paired seeds="<<paired.size()<<"/"<<expected_seeds<<".</p><table><tr><th>Method</th><th>Process mean (us)</th><th>Ratio to target</th><th>Success</th><th>Fallback</th><th>Support</th><th>Shortened</th></tr>";
        std::vector<double> heatrow;
        for(auto& method:baseline) heatrow.push_back(mean(method)/tm); heatrow.push_back(best/tm);
        heat[m[0]+(erasure?"_known_erasures":"_unknown_errors")].push_back({m[14],heatrow});
        for(auto& entry:g.methods) {
            auto& name=entry.first; auto& t=entry.second; bool ablation=name.rfind("hybrid_",0)==0;
            methods<<std::setprecision(17)<<m[0]<<','<<m[14]<<','<<name<<','<<t.mean()<<','<<quantile(t.us,.5)<<','<<quantile(t.us,.95)<<','<<t.process.size()<<','<<t.distinct.size()<<','<<t.calls<<','<<t.fraction(t.correct)<<','<<t.fraction(t.fall)<<','<<t.fraction(t.support)<<','<<t.fraction(t.shortened)<<','<<t.fraction(t.zero_branch)<<','<<t.fraction(t.local)<<','<<t.fraction(t.first_good)<<','<<t.fraction(t.first_captured)<<','<<t.fraction(t.first_both)<<','<<ablation<<'\n';
            std::string color="missing"; double ratio=NAN;
            bool has_paired=true; for(auto s:paired) has_paired&=t.process.count(s);
            if(has_paired&&t.correct==t.calls) {
                ratio=mean(name)/tm; color=ratio>1?"green":"gray";
                if(name!=target) {
                    auto single_ci=bootstrap(g,target,{name},paired,resamples,stable_hash(m[0]+m[14]+name));
                    ratios<<m[0]<<','<<m[14]<<','<<target<<','<<name<<','<<mean(name)<<','<<tm<<','<<ratio<<',';
                    if(std::isfinite(single_ci.first)) ratios<<single_ci.first; ratios<<',';
                    if(std::isfinite(single_ci.second)) ratios<<single_ci.second;
                    ratios<<','<<color<<','<<paired.size()<<','<<ablation<<'\n';
                }
            }
            html<<"<tr><td>"<<escaped(name)<<"</td><td>"<<t.mean()<<"</td><td class='"<<color<<"'>";
            if(std::isfinite(ratio)) html<<ratio; else html<<"n/a";
            html<<"</td><td>"<<t.fraction(t.correct)<<"</td><td>"<<t.fraction(t.fall)<<"</td><td>"<<t.fraction(t.support)<<"</td><td>"<<t.fraction(t.shortened)<<"</td></tr>";
        }
        html<<"<tr><th>Minimum of required baselines</th><td>"<<best<<"</td><td class='"<<(best>tm?"green":"gray")<<"'>"<<best/tm<<"</td><td colspan='4'>"<<escaped(best_name)<<"; 95% paired seed interval: "<<ci.first<<" to "<<ci.second<<"</td></tr></table>";
    }
    std::filesystem::create_directories(outdir/"figures");
    for(auto& entry:heat) {
        long h=135+24*static_cast<long>(entry.second.size()); std::ofstream svg(outdir/"figures"/(entry.first+".svg"));
        svg<<"<svg xmlns='http://www.w3.org/2000/svg' width='1300' height='"<<h<<"' viewBox='0 0 1300 "<<h<<"'><rect width='100%' height='100%' fill='white'/><g font-family='sans-serif' font-size='12'><text x='12' y='22'>"<<escaped(entry.first)<<" : comparator mean / hybrid mean</text><text x='12' y='42'>Green &gt; 1; gray &lt;= 1. Known erasure rows use global/local recovery. Colors do not imply significance.</text>";
        std::vector<std::string> names{"GF","GS","RA","LA","IS","PI","MR","Gao C++","minimum"};
        for(size_t j=0;j<names.size();++j) svg<<"<text x='"<<290+j*108<<"' y='76'>"<<names[j]<<"</text>";
        long row=0;
        for(auto& line:entry.second) {
            long y=90+24*row++; svg<<"<text x='12' y='"<<y+16<<"'>"<<escaped(line.first)<<"</text>";
            for(size_t j=0;j<line.second.size();++j) {
                // Erasure rows have only global/local and its minimum; place in
                // first and last columns, not under an unrelated algorithm.
                size_t col=line.second.size()==2?(j==0?0:8):j;
                auto ratio=line.second[j]; svg<<"<rect x='"<<280+col*108<<"' y='"<<y<<"' width='104' height='22' fill='"<<(ratio>1?"#cdebd1":"#dedede")<<"'/><text x='"<<288+col*108<<"' y='"<<y+16<<"'>"<<std::setprecision(5)<<ratio<<"</text>";
            }
        } svg<<"</g></svg>";
        html<<"<p><a href='figures/"<<escaped(entry.first)<<".svg'>Complete heatmap: "<<escaped(entry.first)<<"</a></p>";
    }
    html<<"<p>Separate code-family controls: each worker's gabidulin_controls.csv. Exact theoretical counts and their probabilities: theory_exact_counts.tsv and theory.csv. Raw sent words, errors and received words are recorded in inputs.txt. No original Sage timings are mixed into these results.</p>";
    std::cout<<"Report written to "<<outdir<<'\n';
}
} // namespace hhy
