#include <cstdio>
#include <memory>
#include <cmath>
#include <vector>
#include <sys/stat.h>

int main(int argc, char **argv) {
    std::FILE *ifp = std::fopen(argv[1], "rb");
    struct stat st;
    ::stat(argv[1], &st);
    const size_t nelem = std::ceil(std::sqrt(st.st_size / 4));
    const size_t npairs = st.st_size / 4;
    std::fprintf(stderr, "Pairwise distances between %zu items\n", nelem);
    size_t nperbucket = 64, nbuckets = 256;
    if(argc > 2) {
        nperbucket = std::strtoull(argv[2], 0, 10);
        if(argc > 3) {
            nbuckets = std::strtoull(argv[3], 0, 10);
        }
    }
    std::fprintf(stderr, "Bucket: %zu/%zu\n", nperbucket, nbuckets);
    const float mul = float(nbuckets);
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> vals;
    vals.resize(nbuckets + 2);
    size_t idx = 0;
    size_t nbuckets_unfilled = vals.size();
    float v;
    for(size_t i = 0; i < nelem; ++i) {
        for(size_t j = i + 1; j < nelem; ++j, ++idx) {
            std::fread(&v, sizeof(float), 1, ifp);
            auto bktid = v == 0.f ? 0: v == 1.f ? int32_t(vals.size() - 1): static_cast<int32_t>(v * mul) + 1;
            auto &vbkt = vals[bktid];
            if(vbkt.size() < nperbucket) {
                vbkt.push_back({i, j});
                if(vbkt.size() == nperbucket) {
                    if(--nbuckets_unfilled == 0) {
                        std::fprintf(stderr, "All buckets are full, early stopping\n");
                        goto end;
                    }
                    std::fprintf(stderr, "[%zu/%zu] %zu buckets left at < %zu entries\n", idx, npairs, nbuckets_unfilled, nperbucket);
                }
            }
        }
    }
    std::fprintf(stderr, "EOF at %zu\n", idx);
    end:
    for(size_t i = 0; i < vals.size(); ++i) {
        if(i == 0)
            std::fprintf(stdout, "Bucket %zu [0.]", i);
        else if(i == vals.size() - 1)
            std::fprintf(stdout, "Bucket %zu [1.]", i);
        else
            std::fprintf(stdout, "Bucket %zu [%g->%g]", i, (i - 1) / double(nbuckets), i / double(nbuckets));
        auto &vbkt = vals[i];
        for(size_t j = 0; j < vbkt.size(); ++j) {
            auto [xi, yi] = vbkt[j];
            std::fprintf(stdout, "\t%u:%u", xi, yi);
        }
        std::fputc('\n', stdout);
    }
    std::fclose(ifp);
}
