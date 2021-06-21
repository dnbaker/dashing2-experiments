#include <cstdio>
#include <memory>
#include <vector>

int main(int argc, char **argv) {
    std::FILE *ifp = std::fopen(argv[1], "rb");
    size_t bufsz = 1 << 18, nelem = bufsz / 4;
    size_t nperbucket = 64, nbuckets = 256;
    if(argc > 2) {
        nperbucket = std::strtoull(argv[2], 0, 10);
        if(argc > 3) {
            nbuckets = std::strtoull(argv[3], 0, 10);
        }
    }
    std::fprintf(stderr, "Bucket: %zu/%zu\n", nperbucket, nbuckets);
    const float mul = float(nbuckets);
    std::unique_ptr<char[]> buf(new char[bufsz]);
    float *ptr = (float *)buf.get();
    std::vector<std::vector<uint64_t>> vals;
    vals.resize(nbuckets + 1);
    size_t idx = 0;
    size_t nbuckets_unfilled = vals.size();
    float v;
    for(;;) {
        if(std::fread(&v, sizeof(float), 1, ifp) != 1) break;
        auto bktid = v == 1.f ? int32_t(vals.size() - 1): static_cast<int32_t>(v * mul);
        auto &vbkt = vals[bktid];
        if(vbkt.size() < nperbucket) {
            vbkt.push_back(idx);
            if(vbkt.size() == nperbucket) {
                if(--nbuckets_unfilled == 0) {
                    std::fprintf(stderr, "All buckets are full, early stopping\n");
                }
                std::fprintf(stderr, "%zu buckets left at < %zu entries\n", nbuckets_unfilled, nperbucket);
            }
        }
        ++idx;
    }
    std::fprintf(stderr, "EOF at %zu\n", idx);
    end:
    for(size_t i = 0; i < vals.size(); ++i) {
        std::fprintf(stdout, "Bucket %zu [%g->%g]", i, i / double(nbuckets), (i + 1) / double(nbuckets));
        auto &vbkt = vals[i];
        for(size_t j = 0; j < vbkt.size(); ++j)
            std::fprintf(stdout, "\t%zu", size_t(vbkt[j]));
        std::fputc('\n', stdout);
    }
    std::fclose(ifp);
}
