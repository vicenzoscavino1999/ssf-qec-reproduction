// Built against the unchanged sector engine and the paper's two seed models.
#define main reference_census_main
#include "sector_reference.cpp"
#undef main

struct Tally {
    U total = 0, flag = 0, input_nonzero = 0, silent = 0;
};

void observe(Engine &g, const int *sites, int k, Tally &counts) {
    g.clear();
    Word readout;
    for (int a = 0; a < k; ++a) {
        const int v = sites[a];
        if (v < g.n) {
            g.flip(g.data[v]);
        } else {
            const int j = v - g.n;
            g.flip(g.measure[j]);
            readout.toggle(j);
        }
    }
    const int start = g.obs.size();
    bool moved = false;
    while (I win = g.winner()) {
        moved = true;
        g.flip(g.acts[g.nacts - (win & g.low)]);
    }
    require(moved == (g.obs.size() < start), "strict descent");
    ++counts.total;
    counts.input_nonzero += (start != 0);
    counts.flag += (moved && g.obs.size() != 0);
    counts.silent += (g.obs == readout && g.tag != 0);
}

int main(int argc, char **argv) {
    if (argc != 3) return 1;
    Engine g(argv[1]);
    const std::string mode(argv[2]);
    if (mode == "flag") {
        std::array<std::array<Tally, 3>, 3> table{};
        const int sites = g.n + g.checks;
        observe(g, nullptr, 0, table[0][0]);
        for (int i = 0; i < sites; ++i) {
            const int a = (i < g.n);
            observe(g, &i, 1, table[a][1-a]);
            for (int j = i+1; j < sites; ++j) {
                const int pair[2] = {i,j};
                const int e = a + (j < g.n);
                observe(g, pair, 2, table[e][2-e]);
            }
        }
        for (int e = 0; e <= 2; ++e) {
            for (int b = 0; b <= 2-e; ++b) {
                const Tally &v = table[e][b];
                std::cout << e << " " << b << " " << v.total << " " << v.flag
                          << " " << v.input_nonzero << " " << v.silent << "\n";
            }
        }
    } else if (mode == "four" && g.n == 225) {
        U total = 0, fails = 0;
        int e[4];
        for (e[0] = 0; e[0] < g.n; ++e[0]) {
            for (e[1] = e[0]+1; e[1] < g.n; ++e[1]) {
                for (e[2] = e[1]+1; e[2] < g.n; ++e[2]) {
                    for (e[3] = e[2]+1; e[3] < g.n; ++e[3]) {
                        ++total;
                        fails += (g.run(e,4,nullptr,0) != 0);
                    }
                }
            }
        }
        std::cout << total << " " << fails << "\n";
    } else {
        return 1;
    }
    return 0;
}
