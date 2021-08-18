#include <atomic>
#include <memory>
#include <algorithm>
#include <vector>
#include <cassert>
#include <cstring>
#include <tuple>
#include <string>
#include <cstdio>
#include <unordered_map>
#include <zlib.h>


template<typename K, typename V>
using map = std::unordered_map<K, V>;

int str2ms(const char *s) {
    int ret =0;
    switch(s[0]) {
        case 'A': case 'a': ret = 0; break;
        case 'C': case 'c': ret = 1; break;
        case 'G': case 'g': ret = 2; break;
        case 'T': case 't': ret = 3; break;
        default: ret = 15;
    }
    ret <<= 4;
    switch(s[1]) {
        case 'A': case 'a': ret |= 0; break;
        case 'C': case 'c': ret |= 1; break;
        case 'G': case 'g': ret |= 2; break;
        case 'T': case 't': ret |= 3; break;
        default: ret |= 15;
    }
    return ret;
}

static inline bool matchchr(const char *s) {
    // | 32 lower-cases the letter, allowing C and c to match
    return (*s | char(32)) == 'c' && s[1] == 'h' && s[2] == 'r';
}

auto parse_file(gzFile ifp) {
    std::vector<uint32_t> contigids;
    std::vector<uint64_t> ids;
    std::vector<uint32_t> starts, stops;
    std::vector<uint8_t> startms, stopms;
    std::atomic<uint64_t> idcounter;
    std::vector<uint16_t> counts;
    std::vector<uint64_t> indptr;
    indptr.push_back(0);
    map<std::string, uint32_t> contignames;
    idcounter.store(0);
    std::string cname;
    size_t ln = 0;
    std::unique_ptr<char[]> buf(new char[1<<20]);
    ssize_t rc;
    for(char *lptr; (lptr = gzgets(ifp, lptr, 1ull << 20)) != nullptr; ++ln) {
        if(ln % 65536 == 0) std::fprintf(stderr, "Processed %zu lines, last rc is %zd\n", ln, rc);
        const uint64_t myid = idcounter++;

        char *p = std::strchr(lptr, '\t');
        assert(p);
        char *p2 = std::strchr(p + 1, '\t');
        assert(p2);
        if(matchchr(p)) p += 3;
        cname.assign(p, p2);
        uint32_t mycid = cname.size();
        auto cnit = contignames.find(cname);
        if(cnit != contignames.end()) mycid = cnit->second;
        else contignames.emplace(cname, mycid);
        contigids.push_back(mycid);
        p = p2 + 1;
        ssize_t srpos = std::strtoll(p, &p, 10);
        ssize_t sppos = std::strtoll(p, &p, 10);
        p = std::strchr(p + 1, '\t'); // Skip the length field
        assert(p[1] == '+' || p[1] == '-' || p[1] == '?' || !std::fprintf(stderr, "p: %s, %d\n", p, p[1]));
        //p += 3;
        p = std::strchr(p + 3, '\t') + 1;
        startms.push_back(str2ms(p));
        p = std::strchr(p, '\t') + 1;
        stopms.push_back(str2ms(p));
        p = std::strchr(std::strchr(std::strchr(p, '\t') + 1, '\t') + 1, '\t') + 1; // Skip two fields
        uint64_t id;
        int32_t ct;
        size_t nids = 0;
        for(char *p3 = p;*p3 && *p3 != '\t';) {
            id = std::strtoull(p3 + 1, &p3, 10);
            ct = std::strtoll(p3 + 1, &p3, 10);
            counts.push_back(ct);
            ids.push_back(id);
            ++nids;
        }
        indptr.push_back(nids + indptr.back());
    }
    return std::make_tuple(contigids, contignames, counts, ids, indptr, ln);
}

int main(int argc, char **argv) {
    gzFile ifp;
    if(argc == 1 || std::strcmp(argv[1], "-") == 0 || std::strcmp(argv[1], "/dev/stdin") == 0) {
        ifp = gzdopen(STDIN_FILENO, "r");
    } else {
        ifp = gzopen(argv[1], "r");
    }
    std::string outpref = "parsed";
    if(argc > 2) outpref = argv[2];
    if(std::find_if(argv, argv + argc, [](auto x) {return !(std::strcmp("-h", x) && std::strcmp("--help", x));}) != argv + argc) {
        std::fprintf(stderr, "recountcsr: This executable parses a tab-delimited, potentially tabix-compressed/indexed, and writes the input to a several files with a prefix {prefix}.cts.u16, {prefix}.ids.u16, {prefix}.indptr.u64, {prefix}.remap");
        std::fprintf(stderr, "This defaults to parsed, but if a second positional argument is provided, it will use it.\n");
        std::fprintf(stderr, "Example (using stdin): `gzip -dc junctions.bgz | recountcsr`\n");
        std::fprintf(stderr, "Example (using zlib): `recountcsr junctions.bgz`\n");
        std::fprintf(stderr, "Example: `recountcsr junctions.bgz jnct`, which uses the 'jnct' as the prefix.\n");
        std::exit(1);
    }
    auto [cids, cnames, counts, ids, indptr, nf] = parse_file(ifp);
    gzclose(ifp);

    std::FILE *fp;
    map<uint64_t, uint16_t> mapper;
    for(const auto id: ids) {
        if(mapper.find(id) == mapper.end()) {
            auto oldsz = mapper.size();
            mapper.emplace(id, oldsz);
        }
    }
    assert(mapper.size() < 65536);
    for(size_t i = 0; i < ids.size(); ++i)
        ids[i] = mapper[ids[i]];

    if((fp = std::fopen((outpref + ".cts.u16").data(), "wb")) == nullptr) {
        std::fprintf(stderr, "Failed to write counts to disk...\n");
        std::exit(1);
    }
    std::fwrite(counts.data(), 2, counts.size(), fp);
    std::fclose(fp);
    if((fp = std::fopen((outpref + ".ids.u16").data(), "wb")) == nullptr) {
        std::fprintf(stderr, "Failed to write IDs to disk...\n");
        std::exit(1);
    }
    for(const auto id: ids) {
        uint16_t sid(id);
        std::fwrite(&sid, 2, 1, fp);
    }
    std::fclose(fp);
    fp = std::fopen((outpref + ".indptr.u64").data(), "w");
    std::fwrite(indptr.data(), sizeof(indptr[0]), indptr.size(), fp);
    std::fclose(fp);

    fp = std::fopen((outpref + ".remap").data(), "w");
    for(const auto &pair: mapper)
        std::fprintf(fp, "%zu:%zu\n", size_t(pair.first), size_t(pair.second));
    std::fclose(fp);
    std::fprintf(stderr, "Shape: %zu, %zu;%zu nnz\n", mapper.size(), nf, counts.size());

    return 0;
}
