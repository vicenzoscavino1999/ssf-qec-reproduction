// Conditional Bernoulli sampling and census, both with the published native SSF.
#define main splitting_225_original_main
#include "../SPLITTING_HGP_225_VV_20261124/motor_225.cpp"
#undef main

namespace {
constexpr int ndata = 225, nread = 108, draws = 1000000;
constexpr U64 sample_seed = 202611242151ULL;

// Rejection makes the bounded draws uniform and independent of the C++ library's
// implementation of uniform_int_distribution. The zero threshold handles powers of two.
int uniform_index(std::mt19937_64& engine, int size) {
    const U64 bound = static_cast<U64>(size);
    const U64 threshold = (U64(0) - bound) % bound;
    U64 value;
    do { value = engine(); } while (value < threshold);
    return static_cast<int>(value % bound);
}

struct BernoulliSites {
    std::array<Outcome, ndata> data{};
    std::array<Outcome, nread> readout{};

    explicit BernoulliSites(const Geometry& geometry) {
        for (int i = 0; i < ndata; ++i) {
            int found = 0;
            for (const auto& label : geometry.outcomes[i]) {
                if (label.offset.weight() == 0 && label.obs.weight() != 0) {
                    data[i] = label;
                    ++found;
                }
            }
            const int degree = data[i].obs.weight();
            check(found == 1 && (degree == 3 || degree == 4), "clean data label");
        }
        for (int j = 0; j < nread; ++j) {
            const auto& labels = geometry.outcomes[ndata + j];
            check(labels.size() == 2, "readout alphabet");
            readout[j] = labels[1];
            const Word bit = j < 64 ? Word{U64(1) << j, 0} : Word{0, U64(1) << (j-64)};
            check(readout[j].obs == bit && readout[j].offset == bit && readout[j].tag == 0,
                  "readout label");
        }
    }
};

int sample(Geometry& geometry, const BernoulliSites& siteset) {
    std::mt19937_64 engine(sample_seed);
    int hits = 0;
    for (int i = 0; i < draws; ++i) {
        const int first = uniform_index(engine, ndata);
        int second = uniform_index(engine, ndata-1);
        if (second >= first) ++second;
        const int check_site = uniform_index(engine, nread);
        hits += geometry.silent(siteset.data[first] ^ siteset.data[second] ^
                                siteset.readout[check_site]);
    }
    return hits;
}

int census(Geometry& geometry, const BernoulliSites& siteset) {
    int hits = 0;
    for (int first = 0; first < ndata; ++first) {
        for (int second = first+1; second < ndata; ++second) {
            const auto pair = siteset.data[first] ^ siteset.data[second];
            for (const auto& reading : siteset.readout)
                hits += geometry.silent(pair ^ reading);
        }
    }
    return hits;
}
}

int main(int argc, char** argv) {
    try {
        check(argc == 2, "use: estrato_225 sample|census");
        const std::string mode = argv[1];
        check(mode == "sample" || mode == "census", "unrecognized mode");
        const auto start = Clock::now();
        Geometry geometry;
        const BernoulliSites siteset(geometry);
        const auto loop_start = Clock::now();
        const double load_seconds = std::chrono::duration<double>(loop_start-start).count();
        const int hits = mode == "sample" ? sample(geometry, siteset) : census(geometry, siteset);
        const double loop_seconds = seconds(loop_start);
        const U64 expected = mode == "sample" ? draws : ndata*(ndata-1)/2*nread;
        check(geometry.queries == expected, "wrong oracle call count");
        std::cout << std::setprecision(12)
                  << "{\"mode\":\"" << mode << "\",\"hits\":" << hits
                  << ",\"calls\":" << geometry.queries << ",\"winner_calls\":" << geometry.winner_calls
                  << ",\"seed\":" << (mode == "sample" ? std::to_string(sample_seed) : "null")
                  << ",\"load_s\":" << load_seconds << ",\"evaluation_s\":" << loop_seconds << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
