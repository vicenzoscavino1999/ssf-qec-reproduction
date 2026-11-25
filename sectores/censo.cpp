#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib")
#else
#include <sys/resource.h>
#endif

using U = std::uint64_t;
using I = std::uint32_t;
using H = std::uint16_t;
using Clock = std::chrono::steady_clock;
double since(Clock::time_point t) {
    return std::chrono::duration<double>(Clock::now() - t).count();
}
void require(bool ok, const char *msg) {
    if (!ok)
        throw std::runtime_error(msg);
}
U peak() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS m{};
    require(GetProcessMemoryInfo(GetCurrentProcess(), &m, sizeof(m)) != 0, "memory");
    return m.PeakWorkingSetSize;
#else
    rusage m{};
    getrusage(RUSAGE_SELF, &m);
#if defined(__APPLE__)
    return m.ru_maxrss;
#else
    return 1024ULL * m.ru_maxrss;
#endif
#endif
}
struct Input {
    std::ifstream f;
    explicit Input(const std::string &file) : f(file, std::ios::binary) {
        require(bool(f), "model file");
    }
    template <class T> T read() {
        T q{};
        f.read(reinterpret_cast<char *>(&q), sizeof(q));
        require(bool(f), "model read");
        return q;
    }
    std::vector<H> cells() {
        int n = read<std::uint8_t>();
        std::vector<H> out(n);
        for (auto &i : out)
            i = read<H>();
        return out;
    }
};
struct Word {
    std::array<U, 7> v{};
    void toggle(int b) { v[b >> 6] ^= U(1) << (b & 63); }
    int size() const {
        int n = 0;
        for (U a : v) {
#ifdef _MSC_VER
            n += int(__popcnt64(a));
#else
            n += __builtin_popcountll(a);
#endif
        }
        return n;
    }
    bool operator==(const Word &o) const { return v == o.v; }
};
struct Shape {
    Word w;
    U tag = 0;
    std::vector<std::pair<H, H>> effects;
};
struct Engine {
    int n, checks, gens;
    I nacts, low;
    std::vector<std::array<H, 12>> rects;
    std::vector<std::vector<std::pair<H, H>>> readers;
    std::vector<I> tbl;
    std::vector<Shape> acts, data, measure;
    std::array<std::vector<unsigned char>, 4096> choices;
    std::vector<I> byshape;
    std::vector<unsigned char> weights;
    int nshapes = 0;
    std::vector<std::pair<int, int>> pairs;
    std::array<H, 432> local{}, used{};
    std::array<U, 432> stamp{};
    int nu = 0;
    U generation = 0, calls = 0, competitions = 0, probes = 0;
    Word obs;
    U tag = 0;
    Shape shape(const std::vector<H> &cells, U mark) {
        Shape a;
        a.tag = mark;
        std::map<H, H> changes;
        for (H b : cells) {
            a.w.toggle(b);
            for (auto [g, bit] : readers[b])
                changes[g] ^= bit;
        }
        for (auto [g, value] : changes)
            if (value)
                a.effects.emplace_back(g, value);
        return a;
    }
    void clear() {
        ++generation;
        nu = 0;
        obs = Word{};
        tag = 0;
    }
    void flip(const Shape &s) {
        for (int i = 0; i < 7; ++i)
            obs.v[i] ^= s.w.v[i];
        tag ^= s.tag;
        for (auto [g, b] : s.effects) {
            if (stamp[g] != generation) {
                stamp[g] = generation;
                used[nu++] = g;
                local[g] = 0;
            }
            local[g] ^= b;
        }
    }
    I winner() {
        ++competitions;
        probes += nu;
        I result = 0;
        for (int i = 0; i < nu; ++i) {
            const I g = used[i];
            result = std::max(result, tbl[4096 * g + local[g]]);
        }
        return result;
    }
    U run(const int *e, int a, const int *d, int b) {
        ++calls;
        clear();
        for (int j = 0; j < b; ++j)
            flip(measure[d[j]]);
        Word goal = obs;
        if (b >= 2 && winner())
            return 0;
        for (int j = 0; j < a; ++j)
            flip(data[e[j]]);
        while (obs.size() > b) {
            I best = winner();
            if (!best)
                return 0;
            flip(acts[nacts - (best & low)]);
        }
        return obs == goal ? tag : 0;
    }
    explicit Engine(const std::string &path) {
        Input in(path);
        std::string id;
        for (int i = 0; i < 8; ++i)
            id += in.read<char>();
        require(id == "SSFMARK1", "magic");
        n = in.read<I>();
        checks = in.read<I>();
        gens = in.read<I>();
        nacts = in.read<I>();
        low = (I(1) << in.read<I>()) - 1;
        int np = in.read<I>();
        require(checks <= 448 && gens <= 432, "dimensions");
        rects.resize(gens);
        readers.resize(checks);
        for (int g = 0; g < gens; ++g)
            for (int j = 0; j < 12; ++j) {
                H b = in.read<H>();
                rects[g][j] = b;
                readers[b].emplace_back(H(g), H(1 << j));
            }
        tbl.resize(4096 * gens);
        in.f.read(reinterpret_cast<char *>(tbl.data()), tbl.size() * sizeof(I));
        for (auto *v : {&acts, &data}) {
            unsigned count = v == &acts ? nacts : unsigned(n);
            v->reserve(count);
            for (unsigned i = 0; i < count; ++i) {
                auto cc = in.cells();
                auto tagval = in.read<U>();
                v->push_back(shape(cc, tagval));
            }
        }
        for (int i = 0; i < checks; ++i)
            measure.push_back(shape({H(i)}, 0));
        for (int i = 0; i < np; ++i) {
            auto a = in.read<H>();
            auto b = in.read<H>();
            pairs.emplace_back(a, b);
        }
        int nt = in.read<I>();
        for (int i = 0; i < nt; ++i) {
            auto s = in.cells(), end = in.cells();
            U expected = in.read<U>();
            clear();
            for (int j : s)
                flip(measure[j]);
            while (I win = winner())
                flip(acts[nacts - (win & low)]);
            Word w;
            for (int j : end)
                w.toggle(j);
            require(obs == w && tag == expected, "global-scan comparison");
        }
        nshapes = in.read<I>();
        for (auto &row : choices) {
            int size = in.read<std::uint8_t>();
            while (size--)
                row.push_back(in.read<std::uint8_t>());
        }
        byshape.resize(gens * nshapes);
        in.f.read(reinterpret_cast<char *>(byshape.data()), byshape.size() * sizeof(I));
        weights.resize(nacts);
        in.f.read(reinterpret_cast<char *>(weights.data()), weights.size());
        require(in.f.peek() == std::ifstream::traits_type::eof(), "extra data");
        competitions = probes = 0;
    }
    std::vector<I> tied(I best) {
        std::vector<I> out;
        I rank = best & ~low;
        unsigned weight = weights[nacts - (best & low)];
        for (int i = 0; i < nu; ++i) {
            int g = used[i];
            H v = local[g];
            I k = tbl[4096 * g + v];
            if ((k & ~low) != rank || weights[nacts - (k & low)] != weight)
                continue;
            for (unsigned char j : choices[v]) {
                I ix = byshape[g * nshapes + j];
                if (weights[ix] == weight && std::find(out.begin(), out.end(), ix) == out.end())
                    out.push_back(ix);
            }
        }
        return out;
    }
    bool possible(const Word &goal) {
        if (obs.size() <= 1)
            return obs == goal && tag != 0;
        I win = winner();
        if (!win)
            return false;
        auto actions = tied(win);
        for (I i : actions) {
            flip(acts[i]);
            bool ok = possible(goal);
            flip(acts[i]);
            if (ok)
                return true;
        }
        return false;
    }
    // The union bound separately charges triples whose first move removes
    // true originals to a failing physical pair with this same readout.
    bool exceptional(const int *e, int d) {
        ++calls;
        clear();
        flip(measure[d]);
        Word goal = obs;
        for (int j = 0; j < 3; ++j)
            flip(data[e[j]]);
        if (obs.size() <= 1)
            return obs == goal && tag != 0;
        I win = winner();
        if (!win)
            return false;
        auto all = tied(win);
        for (I i : all) {
            const auto &act = acts[i];
            bool good = false;
            for (int j = 0; j < 3; ++j)
                if (act.tag == data[e[j]].tag && act.w == data[e[j]].w)
                    good = true;
            if (!good) {
                for (int sub : {3, 5, 6, 7}) {
                    Word w;
                    U z = 0;
                    for (int j = 0; j < 3; ++j)
                        if (sub & (1 << j)) {
                            const auto &f = data[e[j]];
                            z ^= f.tag;
                            for (int k = 0; k < 7; ++k)
                                w.v[k] ^= f.w.v[k];
                        }
                    if (w == act.w && z == act.tag) {
                        good = true;
                        break;
                    }
                }
            }
            if (!good) {
                flip(act);
                bool ok = possible(goal);
                flip(act);
                if (ok)
                    return true;
            }
        }
        return false;
    }
};

struct OrbitRow {
    std::array<unsigned char, 4> v;
    I size;
    std::vector<I> perm;
    std::vector<std::array<unsigned char, 3>> members;
};
void orbit_sum(Engine &g, const std::string &file) {
    auto start = Clock::now();
    Input in(file);
    std::string id;
    for (int i = 0; i < 8; ++i)
        id += in.read<char>();
    bool pure = id == "SSFZERO1";
    require(pure || id == "SSFORBT1", "orbit input");
    int A = g.n == 225 ? 12 : 24, B = A * 3 / 4;
    U total = 0, exceptional = 0, tested = 0;
    int e[3];
    std::map<U, U> tags;
    for (int v = 0; v <= 3; ++v) {
        auto t = Clock::now();
        int np = in.read<I>();
        std::array<std::vector<OrbitRow>, 2> axes;
        for (auto &rows : axes) {
            rows.resize(in.read<I>());
            for (auto &row : rows) {
                for (auto &q : row.v)
                    q = in.read<std::uint8_t>();
                row.size = in.read<I>();
                row.perm.resize(np);
                for (auto &x : row.perm)
                    x = in.read<I>();
                if (pure) {
                    row.members.resize(row.size);
                    for (auto &m : row.members)
                        for (auto &p : m)
                            p = in.read<std::uint8_t>();
                }
            }
        }
        for (I i = 0; i < axes[0].size(); ++i) {
            const auto &x = axes[0][i];
            if (*std::min_element(x.perm.begin(), x.perm.end()) < i)
                continue;
            for (I j = 0; j < axes[1].size(); ++j) {
                const auto &y = axes[1][j];
                int st = 0;
                bool keep = true;
                for (int k = 0; k < np; ++k)
                    if (x.perm[k] == i) {
                        if (y.perm[k] < j) {
                            keep = false;
                            break;
                        }
                        if (y.perm[k] == j)
                            ++st;
                    }
                if (!keep)
                    continue;
                for (int k = 0; k < 3; ++k)
                    e[k] = k < v ? A * x.v[k] + y.v[k] : A * A + B * x.v[k] + y.v[k];
                if (e[0] == e[1] || e[0] == e[2] || e[1] == e[2])
                    continue;
                U ordered = U(x.size) * y.size;
                require(st && ordered % st == 0, "unordered orbit mass");
                U mult = ordered / st;
                total += mult;
                ++tested;
                if (!pure) {
                    if (g.exceptional(e, B * x.v[3] + y.v[3]))
                        exceptional += mult;
                    continue;
                }
                g.clear();
                for (int h : e)
                    g.flip(g.data[h]);
                if (!g.possible(Word{}))
                    continue;
                std::map<U, U> local_hits;
                for (const auto &u : x.members)
                    for (const auto &w : y.members) {
                        for (int k = 0; k < 3; ++k)
                            e[k] = k < v ? A * u[k] + w[k] : A * A + B * u[k] + w[k];
                        U mark = g.run(e, 3, nullptr, 0);
                        if (mark)
                            ++local_hits[mark];
                    }
                for (auto [mark, count] : local_hits) {
                    require(count % st == 0, "logical multiplicity");
                    tags[mark] += count / st;
                    exceptional += count / st;
                }
            }
        }
        std::cout << "# vv " << v << " seconds " << since(t) << " representatives " << tested
                  << std::endl;
    }
    require(total == U(g.n) * (g.n - 1) * (g.n - 2) / 6 * (pure ? 1 : g.checks),
            "complete orbit mass");
    std::cout << (pure ? "triples " : "orbits ") << total << " " << exceptional << " " << tested
              << " " << g.competitions << " " << g.probes << " " << since(start) << " " << peak()
              << " " << g.calls << " ";
    for (auto [mark, count] : tags)
        std::cout << mark << ":" << count << ",";
    std::cout << std::endl;
}
int main(int argc, char **argv) {
    try {
        require(argc == 2 || argc == 3,
                "usage: censo model.bin [orbits.bin | --direct | --allpairs1]");
        auto start = Clock::now();
        Engine g(argv[1]);
        std::cout << "# prep_s " << since(start) << " actions " << g.nacts << " pairs "
                  << g.pairs.size() << std::endl;
        std::string mode = argc == 3 ? argv[2] : "";
        if (!mode.empty() && mode != "--direct" && mode != "--allpairs1") {
            orbit_sum(g, mode);
            return 0;
        }
        int A = g.n == 225 ? 12 : 24, B = A * 3 / 4, m = g.checks;
        std::vector<std::vector<std::pair<int, int>>> reads(A);
        if (mode.empty())
            for (int row = 0; row < A; ++row) {
                std::vector<unsigned char> together(m * m);
                for (const auto &rect : g.rects) {
                    bool contains = false;
                    for (H c : rect)
                        if (c / B == row)
                            contains = true;
                    if (contains)
                        for (H c : rect)
                            for (H k : rect)
                                together[c * m + k] = 1;
                }
                for (int c = 0; c < m; ++c)
                    for (int k = 0; k < c; ++k)
                        if (c / B == row || k / B == row || together[c * m + k])
                            reads[row].emplace_back(c, k);
            }
        int e[3]{}, d[2]{};
        for (auto [a, b] : std::array<std::pair<int, int>, 4>{{{2, 0}, {2, 1}, {2, 2}, {3, 0}}}) {
            if (a == 3 && mode != "--direct")
                continue;
            if (b == 2 && mode == "--allpairs1")
                continue;
            auto t = Clock::now();
            auto old = g.calls, w0 = g.competitions, r0 = g.probes;
            std::map<U, U> hits;
            auto visit = [&]() {
                U tag = g.run(e, a, d, b);
                if (tag)
                    ++hits[tag];
            };
            if (a == 3) {
                for (e[0] = 0; e[0] < g.n; ++e[0])
                    for (e[1] = 0; e[1] < e[0]; ++e[1])
                        for (e[2] = 0; e[2] < e[1]; ++e[2])
                            visit();
            } else if (mode == "--allpairs1") {
                for (e[0] = 0; e[0] < g.n; ++e[0])
                    for (e[1] = 0; e[1] < e[0]; ++e[1])
                        if (!b)
                            visit();
                        else
                            for (d[0] = 0; d[0] < m; ++d[0])
                                visit();
            } else
                for (auto [i, j] : g.pairs) {
                    e[0] = i;
                    e[1] = j;
                    int row = i / A;
                    if (!b)
                        visit();
                    else if (b == 1) {
                        int begin = mode.empty() ? B * row : 0, end = mode.empty() ? begin + B : m;
                        for (d[0] = begin; d[0] < end; ++d[0])
                            visit();
                    } else if (mode.empty()) {
                        for (auto [c, k] : reads[row]) {
                            d[0] = c;
                            d[1] = k;
                            visit();
                        }
                    } else
                        for (d[0] = 0; d[0] < m; ++d[0])
                            for (d[1] = 0; d[1] < d[0]; ++d[1])
                                visit();
                }
            U count = 0;
            for (auto [tag, k] : hits)
                count += k;
            if (!b && a == 2)
                require(count == 0, "empty-readout pair assumption");
            std::cout << a << " " << b << " exact " << g.calls - old << " " << count << " "
                      << hits.size() << " " << g.competitions - w0 << " " << g.probes - r0 << " "
                      << since(t) << " ";
            for (auto [mark, k] : hits)
                std::cout << mark << ":" << k << ",";
            std::cout << std::endl;
        }
        std::cout << "# total_s " << since(start) << " peak_bytes " << peak() << std::endl;
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
