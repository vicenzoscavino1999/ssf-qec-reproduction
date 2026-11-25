// Independent draws and conditional Metropolis chains for the affine 225-qubit code.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <immintrin.h>
#endif

using U64 = std::uint64_t;
using U32 = std::uint32_t;
using Clock = std::chrono::steady_clock;
constexpr int sites = 333, generators = 108, seeds = 12, repeats = 3;
constexpr int first_length = 100000, max_length = 3100000;
constexpr double grid = 100000000.0;
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now()-start).count();
}
struct Word {
    U64 lo=0, hi=0;
    bool operator==(Word v) const { return lo==v.lo && hi==v.hi; }
    Word operator^(Word v) const { return {lo^v.lo,hi^v.hi}; }
    int weight() const {
#if defined(_MSC_VER)
        return int(__popcnt64(lo)+__popcnt64(hi));
#else
        return __builtin_popcountll(lo)+__builtin_popcountll(hi);
#endif
    }
};
struct Outcome {
    Word obs, offset;
    U32 tag=0;
    Outcome operator^(const Outcome& v) const {
        return {obs^v.obs, offset^v.offset, tag^v.tag};
    }
};
struct Input {
    std::ifstream file;
    Input():file("modelo.bin",std::ios::binary) { check(bool(file),"abrir modelo"); }
    template <typename T> T get() {
        T v{}; file.read(reinterpret_cast<char*>(&v),sizeof(v));
        check(bool(file),"archivo de modelo incompleto"); return v;
    }
    Word word() { return {get<U64>(),get<U64>()}; }
};
struct Rng {
    std::mt19937_64 gen;
    explicit Rng(U64 seed):gen(seed) {}
    double unit() { return static_cast<double>(gen() >> 11)*(1.0/9007199254740992.0); }
    int index(int n) { return std::uniform_int_distribution<int>(0,n-1)(gen); }
};
struct Geometry {
    struct Rect { Word mask; int shift; };
    std::array<Rect,generators> rect;
    std::vector<U32> lookup = std::vector<U32>(generators*4096);
    std::vector<Word> footprint;
    std::vector<U32> action_tag;
    U32 low_mask;
    std::array<std::vector<Outcome>,sites> outcomes;
    std::vector<std::array<double,2>> rates;
    using Laws = std::array<std::vector<double>,sites>;
    std::vector<Laws> laws;
    U64 queries=0, winner_calls=0;

    unsigned index(int g,Word s) const {
        const auto& r=rect[g];
        return unsigned(_pext_u64(s.lo,r.mask.lo) | (_pext_u64(s.hi,r.mask.hi)<<r.shift));
    }
    U32 winner(Word s) {
        ++winner_calls;
        U32 best=0;
        for(int g=0;g<generators;++g) best=std::max(best,lookup[4096*g+index(g,s)]);
        return best;
    }
    bool terminal(Word d) const {
        if(d.weight()<=1) return true; // The smallest nonzero footprint has size three.
        for(int g=0;g<generators;++g) if(lookup[4096*g+index(g,d)]) return false;
        return true;
    }
    std::pair<Word,U32> decode(Word obs) {
        U32 tag=0;
        while(U32 win=winner(obs)) {
            auto ix=footprint.size()-(win&low_mask);
            obs=obs^footprint[ix]; tag^=action_tag[ix];
        }
        return {obs,tag};
    }
    bool silent(Outcome x) {
        ++queries;
        const int need=x.offset.weight();
        if(x.obs==x.offset) return x.tag && terminal(x.offset);
        if(x.obs.weight()<=need || !terminal(x.offset)) return false;
        while(x.obs.weight()>need) {
            U32 win=winner(x.obs);
            if(!win) return false;
            auto ix=footprint.size()-(win&low_mask);
            x.obs=x.obs^footprint[ix]; x.tag^=action_tag[ix];
        }
        return x.obs==x.offset && x.tag!=0;
    }
    Geometry() {
        Input in;
        char id[8]; in.file.read(id,8);
        check(std::string(id,8)=="SSF22501","formato del modelo");
        U32 count=in.get<U32>(); low_mask=in.get<U32>();
        for(auto& r:rect) { r.mask=in.word(); r.shift=Word{r.mask.lo,0}.weight(); }
        in.file.read(reinterpret_cast<char*>(lookup.data()),lookup.size()*sizeof(U32));
        for(U32 i=0;i<count;++i) {footprint.push_back(in.word());action_tag.push_back(in.get<U32>());}
        U32 test_count=in.get<U32>();
        for(U32 i=0;i<test_count;++i) {
            auto s=in.word(), d=in.word(); U32 tag=in.get<U32>();
            auto result=decode(s);
            check(result.first==d && result.second==tag,"decoder native/Python");
            check(!silent({s,d,tag}) && silent({s,d,tag^1U}),"evento del residuo");
            Word other = d ^ Word{1ULL<<unsigned(i%64),0};
            check(!silent({s,other,tag^1U}),"sindrome fisico incorrecto");
        }
        U32 number=in.get<U32>(); rates.resize(number); laws.resize(number);
        for(auto& r:rates) r={in.get<double>(),in.get<double>()};
        for(auto& row:outcomes) {
            row.resize(in.get<U32>());
            for(auto& v:row) v={in.word(),in.word(),in.get<U32>()};
        }
        for(auto& lev:laws) for(int b=0;b<sites;++b) {
            double norm=0;
            for(size_t k=0;k<outcomes[b].size();++k) {
                double v=in.get<double>(); check(v>0 && v<=1,"probabilidad de fibra");
                lev[b].push_back(v); norm+=v;
            }
            check(std::abs(norm-1)<2e-14,"normalizacion");
        }
        check(in.file.peek()==std::ifstream::traits_type::eof(),"datos adicionales");
        queries=winner_calls=0;
    }
};
using Configuration=std::array<unsigned char,sites>;
Outcome reduce(const Geometry& g,const Configuration& cfg) {
    Outcome x{};
    for(int s=0;s<sites;++s) if(cfg[s]) x=x^g.outcomes[s][cfg[s]];
    return x;
}
// Independent site draws. A common Bernoulli envelope permits skipping zero wires.
struct AnchorSampler {
    Geometry& g;
    double pmax=0, inv_log;
    std::array<double,sites> occupied{};
    std::array<std::vector<double>,sites> cdf;
    explicit AnchorSampler(Geometry& g):g(g) {
        for(int s=0;s<sites;++s) {
            const auto& row=g.laws[0][s];
            for(size_t k=1;k<row.size();++k) { occupied[s]+=row[k]; cdf[s].push_back(occupied[s]); }
            for(double& x:cdf[s]) x/=occupied[s];
            cdf[s].back()=1; pmax=std::max(pmax,occupied[s]);
        }
        inv_log=1/std::log1p(-pmax);
    }
    bool draw(Rng& rng, Configuration* saved) {
        Outcome x{};
        std::array<std::pair<int,int>,sites> active;
        int used=0, p=-1;
        for(;;) {
            double gap=std::log1p(-rng.unit())*inv_log;
            if(gap>=sites) break;
            p+=1+int(gap); if(p>=sites) break;
            if(rng.unit()*pmax>=occupied[p]) continue;
            const auto& row=cdf[p];
            int k=1+int(std::upper_bound(row.begin(),row.end(),rng.unit())-row.begin());
            x=x^g.outcomes[p][k];
            if(saved) active[used++]={p,k};
        }
        if(!g.silent(x)) return false;
        if(saved) {
            saved->fill(0);
            for(int i=0;i<used;++i) (*saved)[active[i].first]=static_cast<unsigned char>(active[i].second);
        }
        return true;
    }
};
void anchor(Geometry& g) {
    auto start=Clock::now();
    AnchorSampler sampler(g);
    constexpr int calibration=250000, desired=1600;
    Rng pilot(202611241107ULL), rng(202611241137ULL);
    int pilot_hits=0;
    for(int i=0;i<calibration;++i) pilot_hits+=sampler.draw(pilot,nullptr);
    U64 planned=std::clamp<U64>((U64(desired)*(calibration+1)+pilot_hits)/(pilot_hits+1),1000000,200000000);
    std::cout<<"Calibracion independiente: "<<pilot_hits<<"/"<<calibration
             <<", ancla fija "<<planned<<" historias."<<std::endl;
    std::vector<Configuration> starts;
    int hits=0;
    U64 after_pilot=g.queries;
    for(U64 i=1;i<=planned;++i) {
        Configuration cfg;
        bool save=starts.size()<seeds;
        if(sampler.draw(rng,save?&cfg:nullptr)) {
            ++hits;
            if(save) {check(g.silent(reduce(g,cfg)),"semilla del ancla"); starts.push_back(cfg);}
        }
        if(i%5000000==0) std::cout<<"Ancla "<<i<<"/"<<planned<<", "<<hits<<" fallos."<<std::endl;
    }
    std::ofstream out("ancla.txt"); out<<std::setprecision(17);
    out<<calibration<<" "<<pilot_hits<<" "<<planned<<" "<<hits<<" "<<seconds(start)<<" "
       <<after_pilot<<" "<<g.queries<<" "<<g.winner_calls<<"\n";
    for(const auto& cfg:starts) {for(auto v:cfg) out<<int(v)<<" ";out<<"\n";}
    check(bool(out),"guardar ancla");
    std::cout<<"Ancla completa: "<<hits<<"/"<<planned<<" en "<<seconds(start)<<" s."<<std::endl;
}
struct Bridge {
    std::array<std::vector<double>,sites> field;
    double base=0;
    Bridge(const Geometry::Laws& high,const Geometry::Laws& low) {
        for(int b=0;b<sites;++b) {
            double zero=std::log(high[b][0]/low[b][0]); base+=zero;
            for(size_t k=0;k<high[b].size();++k) field[b].push_back(std::log(high[b][k]/low[b][k])-zero);
        }
    }
    double score(const Configuration& cfg) const {
        double result=base;
        for(int b=0;b<sites;++b) if(cfg[b]) result+=field[b][cfg[b]];
        return result;
    }
};
struct Change {int start; double minus,plus;};
struct Moments {double sum=0,square=0;};
double logistic(double z) {return 1/(1+std::exp(z));}
struct Summary {
    // [left/right half][previous/next level].
    Moments halves[2][2];
    std::map<long long,int> hist[2];
    double rule=0;
};
Summary summarize(const std::vector<Change>& runs,int n,bool left,bool right,bool bins) {
    Summary s;
    for(size_t i=0;i<runs.size();++i) {
        int a=runs[i].start, b=i+1<runs.size()?runs[i+1].start:n;
        if(a>=b) continue;
        int dwell=b-a;
        int first=std::max(0,std::min(b,n/2)-a);
        int count[2]={first,dwell-first};
        for(int side=0;side<2;++side) {
            if((side==0 && !left) || (side==1 && !right)) continue;
            double z=side?runs[i].plus:runs[i].minus;
            double f=logistic(side?z:-z);
            for(int h=0;h<2;++h) {
                s.halves[h][side].sum+=count[h]*f;
                s.halves[h][side].square+=count[h]*f*f;
            }
            if(bins) s.hist[side][std::llround(z*grid)]+=dwell;
        }
    }
    double deviation=0,sigma=0;
    for(int side=0;side<2;++side) {
        if((!left && !side)||(!right && side)) continue;
        const auto& a=s.halves[0][side]; const auto& b=s.halves[1][side];
        double mean=(a.sum+b.sum)/n;
        double var=std::max(0.0,(a.square+b.square-n*mean*mean)/(n-1));
        sigma=std::max(sigma,std::sqrt(var/n)/mean);
        deviation=std::max(deviation,std::abs(2*a.sum/n-mean)/mean);
    }
    s.rule=sigma+deviation;
    return s;
}
void chains(Geometry& g) {
    std::ifstream init("ancla.txt"); std::string header; std::getline(init,header);
    std::vector<Configuration> starts(seeds);
    for(auto& s:starts) {
        for(auto& v:s) {int c;check(bool(init>>c),"faltan semillas positivas");v=static_cast<unsigned char>(c);}
        check(g.silent(reduce(g,s)),"semilla no silenciosa");
    }
    const int nlev=int(g.laws.size());
    const double tolerance=0.05/std::sqrt(double(nlev));
    std::vector<Bridge> edges;
    for(int j=0;j+1<nlev;++j) edges.emplace_back(g.laws[j],g.laws[j+1]);
    std::ofstream chains("cadenas.tsv"), hist("histogramas.tsv");
    chains<<std::setprecision(17);
    chains<<"seed\trep\tj\tn\tmoves\toracles\trule\tleft0\tleft0sq\tleft1\tleft1sq\tright0\tright0sq\tright1\tright1sq\n";
    hist<<"seed\trep\tj\tside\tscore_1e8\tcount\n";
    auto start=Clock::now();
    for(int l=0;l<seeds;++l) for(int rep=0;rep<repeats;++rep) {
        Configuration cfg=starts[l]; Outcome x=reduce(g,cfg);
        for(int j=0;j<nlev;++j) {
            const bool left=j>0,right=j+1<nlev;
            Rng rng(202611249000ULL+1000000*l+1000*rep+j);
            const auto& law=g.laws[j]; U64 before=g.queries;
            double zm=left?edges[j-1].score(cfg):0, zp=right?edges[j].score(cfg):0;
            std::vector<Change> records{{0,zm,zp}};
            int moves=0,used=0,target=first_length;
            Summary summary;
            for(;;) {
                for(int t=used;t<target;++t) {
                    int s=rng.index(sites),old=cfg[s],k=rng.index(int(law[s].size())-1);
                    if(k>=old) ++k;
                    if(law[s][k]<law[s][old] && rng.unit()*law[s][old]>=law[s][k]) continue;
                    auto trial=x^g.outcomes[s][old]^g.outcomes[s][k];
                    if(!g.silent(trial)) continue;
                    ++moves; cfg[s]=static_cast<unsigned char>(k); x=trial;
                    if(left) zm+=edges[j-1].field[s][k]-edges[j-1].field[s][old];
                    if(right) zp+=edges[j].field[s][k]-edges[j].field[s][old];
                    if(moves%128==0) {
                        if(left) {double v=edges[j-1].score(cfg); check(std::abs(v-zm)<1e-9,"score anterior");zm=v;}
                        if(right) {double v=edges[j].score(cfg);check(std::abs(v-zp)<1e-9,"score siguiente");zp=v;}
                    }
                    // Samples are post-step; a move at t supplies the t-th zero-indexed sample.
                    records.push_back({t,zm,zp});
                }
                used=target; summary=summarize(records,used,left,right,false);
                if(summary.rule<=tolerance || used==max_length) break;
                target=std::min(max_length,2*used+first_length);
            }
            summary=summarize(records,used,left,right,true);
            chains<<l<<'\t'<<rep<<'\t'<<j<<'\t'<<used<<'\t'<<moves<<'\t'<<g.queries-before<<'\t'<<summary.rule;
            for(int side=0;side<2;++side) for(int h=0;h<2;++h) {
                const auto& m=summary.halves[h][side];
                chains<<'\t'<<m.sum<<'\t'<<m.square;
            }
            chains<<'\n';
            for(int side=0;side<2;++side) for(auto [key,n]:summary.hist[side])
                hist<<l<<'\t'<<rep<<'\t'<<j<<'\t'<<side<<'\t'<<key<<'\t'<<n<<'\n';
            check(g.silent(reduce(g,cfg)),"salida del nivel no silenciosa");
        }
        chains.flush();hist.flush();
        if(rep+1==repeats) std::cout<<"Semilla "<<l+1<<"/"<<seeds<<" completa, "<<seconds(start)<<" s."<<std::endl;
    }
    check(bool(chains)&&bool(hist),"salidas de las cadenas");
}
int main(int argc,char** argv) {
    try {
        check(argc==2,"indicar ancla o cadenas");
        Geometry geometry;
        std::string mode=argv[1];
        if(mode=="ancla") anchor(geometry);
        else if(mode=="cadenas") chains(geometry);
        else throw std::runtime_error("modo desconocido");
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
    return 0;
}
