#include <cstdio>
#include <cstring>
#include <memory>
#include <algorithm>
#include <cmath>
#include <vector>
#include <sys/stat.h>

int main(int argc, char **argv) {
    if(std::find_if(argv, argc + argv, [](auto x) {return std::strcmp(x, "-h") == 0 || std::strcmp(x, "--help") == 0;}) != argc + argv) {
        std::fprintf(stderr, "Usage: %s distance_matrix.f32 <nperbucket=64> <nbuckets=256>\n", argv[0]);
        std::exit(1);
    }
    struct stat st;
    if(::stat(argv[1], &st)) {
        std::fprintf(stderr, "input %s is empty.\n", argv[1]);
        std::exit(1);
    }
    const size_t nelem = std::ceil(std::sqrt(st.st_size / sizeof(float) * 2));
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
    std::FILE *ifp = std::fopen(argv[1], "rb");
    ssize_t lastunfilled = -1;
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
                    if(nbuckets_unfilled == 1) {
                        ssize_t ind;
                        if(lastunfilled > 0) ind = lastunfilled;
                        else ind = std::find_if(vals.begin(), vals.end(), [=](const auto &x) {return x.size() < nperbucket;}) - vals.begin();
                        std::fprintf(stderr, "Bucket %zu/%g-%g, Current size %zu/%zu\n", ind, (ind - 1) / double(nbuckets), ind / double(nbuckets), vals[ind].size(), nperbucket);
                    }
                }
            }
        }
    }
    std::fclose(ifp);
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
}
