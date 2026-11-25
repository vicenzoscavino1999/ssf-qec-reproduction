// Upper bound for silent histories with a first move inequivalent to every nonempty subset of the data faults.
#define main reference_census_main
#include "sector_reference.cpp"
#undef main

struct Signature { Word syndrome; U tag = 0; };

bool may_finish(Engine &g, const Word &goal, int readouts, unsigned &budget) {
    if (g.obs.size() <= readouts) return g.obs == goal && g.tag != 0;
    I best = g.winner();
    if (!best) return false;
    if (!budget) return true; // The orbit is retained if the search is unfinished.
    --budget;
    for (I next : g.tied(best)) {
        g.flip(g.acts[next]);
        bool possible = may_finish(g, goal, readouts, budget);
        g.flip(g.acts[next]);
        if (possible) return true;
    }
    return false;
}

bool first_exception(Engine &g, const int *data, int a, const int *d, int b, unsigned budget) {
    g.clear();
    for (int j=0; j<b; ++j) g.flip(g.measure[d[j]]);
    const Word goal = g.obs;
    if (b >= 2 && g.winner()) return false; // D has to be stable at termination.
    for (int j=0; j<a; ++j) g.flip(g.data[data[j]]);
    if (g.obs.size() <= b) return g.obs == goal && g.tag != 0;
    I best = g.winner();
    if (!best) return false;
    std::array<Signature,16> sums{};
    for (int j=0; j<a; ++j) {
        const auto &f = g.data[data[j]];
        for (int h=0; h<(1<<j); ++h) {
            Signature &dest = sums[h+(1<<j)];
            dest.tag = sums[h].tag ^ f.tag;
            for (int k=0; k<7; ++k) dest.syndrome.v[k] = sums[h].syndrome.v[k] ^ f.w.v[k];
        }
    }
    for (I next : g.tied(best)) {
        const Shape &f = g.acts[next];
        bool genuine = false;
        for (int h=1; h<(1<<a); ++h) {
            if (sums[h].tag == f.tag && sums[h].syndrome == f.w) { genuine = true; break; }
        }
        if (genuine) continue; // Charge separately to a smaller residual with the same D.
        g.flip(f);
        bool possible = may_finish(g, goal, b, budget);
        g.flip(f);
        if (possible) return true;
    }
    return false;
}

struct AxisRow { std::array<unsigned char,5> coords{}; I size; std::vector<I> relabel; };
U choose(U n, unsigned k) {
    U out=1;
    for (unsigned j=1; j<=k; ++j) out=out*(n-j+1)/j;
    return out;
}

int main(int argc, char **argv) {
    try {
        require(argc==4,"expected model, orbit file, search budget");
        Engine g(argv[1]); Input in(argv[2]);
        std::string id; for (int k=0; k<8; ++k) id+=in.read<char>();
        require(id=="ABORBIT1","input header");
        int t=in.read<I>(), a=in.read<I>(), b=in.read<I>();
        require((t==1 || t==2) && a>=2 && a<=4 && b>=0 && a+b<=5 && g.n==225*t*t,"dimensions");
        unsigned budget=std::stoul(argv[3]);
        int A=12*t,B=9*t;
        U total=0, exceptional=0, reps=0;
        for (int v=0; v<=a; ++v) {
            I np=in.read<I>();
            std::array<std::vector<AxisRow>,2> axes;
            for (auto &rows:axes) {
                rows.resize(in.read<I>());
                for (auto &r:rows) {
                    for (int h=0; h<a+b; ++h) r.coords[h]=in.read<std::uint8_t>();
                    r.size=in.read<I>();
                    r.relabel.resize(np);
                    for (I &u:r.relabel) u=in.read<I>();
                }
            }
            for (I i=0; i<axes[0].size(); ++i) {
                const auto &x=axes[0][i];
                if (*std::min_element(x.relabel.begin(),x.relabel.end())<i) continue;
                std::vector<I> keep;
                for (I k=0; k<np; ++k) if (x.relabel[k]==i) keep.push_back(k);
                for (I j=0; j<axes[1].size(); ++j) {
                    const auto &y=axes[1][j];
                    unsigned stabilizer=0; bool canonical=true;
                    for (I k:keep) {
                        if (y.relabel[k]<j) { canonical=false; break; }
                        stabilizer+=(y.relabel[k]==j);
                    }
                    if (!canonical) continue;
                    int data[4],readout[3]; bool distinct=true;
                    for (int h=0; h<a; ++h) {
                        data[h]=(h<v)?A*x.coords[h]+y.coords[h]:A*A+B*x.coords[h]+y.coords[h];
                        for(int prev=0; prev<h; ++prev) distinct &= data[h]!=data[prev];
                    }
                    for (int h=0; h<b; ++h) {
                        readout[h]=B*x.coords[a+h]+y.coords[a+h];
                        for(int prev=0; prev<h; ++prev) distinct &= readout[h]!=readout[prev];
                    }
                    if (!distinct) continue;
                    U ordered=U(x.size)*y.size;
                    require(stabilizer && ordered%stabilizer==0,"unordered orbit multiplicity");
                    U mass=ordered/stabilizer;
                    total+=mass; ++reps;
                    if (first_exception(g,data,a,readout,b,budget)) exceptional+=mass;
                }
            }
            std::cerr<<"vv "<<v<<" checked "<<reps<<" orbits"<<std::endl;
        }
        require(in.f.peek()==std::ifstream::traits_type::eof(),"trailing input");
        require(total==choose(g.n,a)*choose(g.checks,b),"orbit partition incomplete");
        std::cout<<a<<" "<<b<<" "<<total<<" "<<exceptional<<" "<<reps<<" "<<budget<<std::endl;
        return 0;
    } catch(const std::exception &e) {std::cerr<<e.what()<<std::endl;return 1;}
}
