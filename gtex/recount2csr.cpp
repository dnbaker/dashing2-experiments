#include <atomic>
#include <ctime>
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
#include <mio.hpp>
#include <getopt.h>
#include <omp.h>


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

void buffer_fp(std::FILE *fp) {
    struct stat st;
    ::fstat(::fileno(fp), &st);
    std::setvbuf(fp, nullptr,  _IOFBF, st.st_blksize);
}

auto parse_file(gzFile ifp, std::string outpref, const size_t bufsize=size_t(1<<26)) {
    std::FILE *ipfp = std::fopen((outpref + ".indptr.u64").data(), "w");
    if(!ipfp) throw 1;
    buffer_fp(ipfp);
    uint64_t ipv = 0;
    std::fwrite(&ipv, sizeof(ipv), 1, ipfp);
    std::FILE *ctfp;
    if((ctfp = std::fopen((outpref + ".cts.u32").data(), "wb")) == nullptr) {
        throw 2;
    }
    buffer_fp(ctfp);
    std::FILE *idfp;
    std::string idpath = outpref + ".ids.u64";
    if((idfp = std::fopen(idpath.data(), "wb")) == nullptr) {
        throw 3;
    }
    buffer_fp(idfp);
    std::vector<uint32_t> contigids;
    //std::vector<uint32_t> starts, stops;
    //std::vector<uint8_t> startms, stopms;
    std::atomic<uint64_t> idcounter{0};
    //idcounter.store(0);
    map<std::string, uint32_t> contignames;
    std::string contig_name;
    size_t ln = 0;
    std::unique_ptr<char[]> buf(new char[bufsize]);
    uint32_t maxct = 0;
    for(char *lptr; (lptr = gzgets(ifp, buf.get(), bufsize)) != nullptr; ++ln) {
        if(ln % 65536 == 0) std::fprintf(stderr, "Processed %zu lines\n", ln);
        const uint64_t myid = idcounter++;

        char *p = std::strchr(lptr, '\t');
        assert(p || !std::fprintf(stderr, "Line %zu failed. line: '%s'\n", ln, lptr));
        char *p2 = std::strchr(p + 1, '\t');
        assert(p2);
        if(matchchr(p)) p += 3;
        contig_name.assign(p, p2);
        uint32_t mycid = contignames.size();
        auto cnit = contignames.find(contig_name);
        if(cnit != contignames.end()) mycid = cnit->second;
        else contignames.emplace(contig_name, mycid);
        contigids.push_back(mycid);
        p = p2 + 1;
        ssize_t srpos = std::strtoll(p, &p, 10);
        ssize_t sppos = std::strtoll(p, &p, 10);
        p = std::strchr(p + 1, '\t'); // Skip the length field
        assert(p[1] == '+' || p[1] == '-' || p[1] == '?' || !std::fprintf(stderr, "p: %s, %d\n", p, p[1]));
        //p += 3;
        p = std::strchr(p + 3, '\t') + 1;
        //startms.push_back(str2ms(p));
        p = std::strchr(p, '\t') + 1;
        //stopms.push_back(str2ms(p));
        p = std::strchr(std::strchr(std::strchr(p, '\t') + 1, '\t') + 1, '\t') + 1; // Skip two fields
        uint64_t id;
        uint32_t ct;
        size_t nids = 0;
        for(char *p3 = p;*p3 && *p3 != '\t';) {
            id = std::strtoull(p3 + 1, &p3, 10);
            ct = std::strtoll(p3 + 1, &p3, 10);
            if(__builtin_expect(ct == 0, 0)) {
                std::fprintf(stderr, "Found a count of 0 in line '%s'\n", lptr);
                std::exit(1);
            }
            std::fwrite(&id, sizeof(id), 1, idfp);
            std::fwrite(&ct, sizeof(ct), 1, ctfp);
            ++nids;
            maxct = std::max(ct, maxct);
        }
        ipv += nids;
        std::fwrite(&ipv, sizeof(ipv), 1, ipfp);
    }
    std::fprintf(stderr, "%zu total lines\n", ln);
    std::fclose(ipfp);
    std::fclose(idfp);
    std::fclose(ctfp);
    return std::make_tuple(contigids, contignames, ln, idpath, (outpref + ".indptr.u64"), (outpref + ".cts.u32"), maxct, ipv);
}

int usage(const char *ex) {
    std::fprintf(stderr, "Usage: %s <junctions.bgz> [optional: outprefix, defaults to 'parsed'\n", ex);
    std::fprintf(stderr, "recountcsr: This executable parses a tab-delimited, potentially tabix-compressed/indexed, and writes the input to a several files with a prefix {prefix}.cts.u16, {prefix}.ids.u16, {prefix}.indptr.u64, {prefix}.remap");
    std::fprintf(stderr, "This defaults to parsed, but if a second positional argument is provided, it will use it.\n");
    std::fprintf(stderr, "Example (using stdin): `gzip -dc junctions.bgz | recountcsr`\n");
    std::fprintf(stderr, "Example (using zlib): `recountcsr junctions.bgz`\n");
    std::fprintf(stderr, "Example: `recountcsr junctions.bgz jnct`, which uses the 'jnct' as the prefix.\n");
    return 1;
}
int getnt() {
    int ret = 0;
    #pragma omp parallel
    {
        ret = omp_get_num_threads();
    }
    return ret;
}

int main(int argc, char **argv) {
    if(std::find_if(argv, argv + argc, [](auto x) {return std::strcmp("--help", x) == 0 || std::strcmp("-h", x) == 0;}) !=  argc + argv) {
        return usage(argv[0]);
    }
    gzFile ifp;
    if(argc == 1 || std::strcmp(argv[1], "-") == 0 || std::strcmp(argv[1], "/dev/stdin") == 0) {
        ifp = gzdopen(STDIN_FILENO, "r");
        if(ifp == nullptr) throw std::runtime_error("Failed to open STDIN fileno");
    } else {
        ifp = gzopen(argv[1], "r");
        if(ifp == nullptr) throw std::runtime_error(std::string("Failed to open ") + argv[1]);
    }
    std::time_t result = std::time(NULL);
    std::string outpref = std::string("parsed") + std::to_string(static_cast<size_t>(result));
    if(argc > 2) {
        outpref = argv[optind + 1];
    }
    auto [cids, cnames, nlines, idpath, indptrpath, ctspath, maxct, nnz] = parse_file(ifp, outpref);
    std::fprintf(stderr, "max ct %u\n", maxct);
    gzclose(ifp);


    size_t njnct = 0, idn = 0; // This is set after filling the mapper
    {
        std::atomic<size_t> idgen{0};
        std::FILE *fp;

        int idfd = ::open(idpath.data(), O_RDWR);
        if(idfd < 0) {
            perror("Failed to open idpath for reading and writing");
            throw std::runtime_error("idpath RDWR O_failure");
        }
        mio::mmap_sink idmm(idfd);
        uint64_t *idptr = (uint64_t *)idmm.data();
        idn = idmm.size() / 8;
        assert(idmm.size() % 8 == 0);
        // Re-name ids in-place, but converting to 32-bits
        if((fp = std::fopen((outpref + ".remap").data(), "w")) == nullptr) {
            throw std::runtime_error(std::string("Failed to open remap file ") + outpref + ".remap");
        }
        const int nt = getnt();
        map<uint64_t, uint64_t> mapper;
        if(nt == 1) {
            std::transform(idptr, &idptr[idn], (uint64_t *)idptr, [&mapper,fp](uint64_t x) -> uint64_t {
                auto it = mapper.find(x);
                if(it == mapper.end()) {
                    it = mapper.emplace(x, mapper.size()).first;
                    std::fprintf(fp, "%zu:%zu\n", size_t(it->first), size_t(it->second));
                }
                return it->second;
            });
        } else {
            #pragma omp parallel for
            for(size_t i = 0; i < idn; ++i) {
                auto id = idptr[i];
                auto it = mapper.find(id);
                if(it == mapper.end()) {
                    #pragma omp critical
                    {
                        if((it = mapper.find(id)) == mapper.end()) {
                            it = mapper.emplace(id, mapper.size()).first;
                        }
                    }
                }
            }
            const auto &cmap = mapper;
            #pragma omp parallel for schedule(static, 8192)
            for(size_t i = 0; i < idn; ++i) {
                auto &res = idptr[i];
                auto it = cmap.find(res);
                if(it == cmap.end()) {std::fprintf(stderr, "Map failed!!!!\n"); std::exit(1);}
                res = it->second;
            }
        }
        ::close(idfd);
        njnct = mapper.size();
        std::fclose(fp);
    }

    std::string idop = outpref + ".ids.u";
    int itype = -1;
    if(njnct <= 255u) {
        idop += "8";
        itype = 0;
    } else if(njnct <= 65535u) {
        idop += "16";
        itype = 1;
    } else if(njnct <= 0xFFFFFFFFu) {
        idop += "32";
        itype = 2;
    } else {
        idop += "64";
        itype = 3;
    }
    if(itype != 3) { // If uint32_t, do nothing!
        int idfd = ::open(idpath.data(), O_RDWR);
        if(idfd < 0) {
            perror("Failed to open idpath for reading and writing");
            throw std::runtime_error("idpath RDWR O_failure");
        }
        mio::mmap_sink idmm(idfd);
        uint64_t *idptr = (uint64_t *)idmm.data();
#define CASE_I(iv, IT)\
        case iv: {\
            std::copy(idptr, idptr + idn, (IT *)idptr);\
            ::ftruncate(idfd, idn * sizeof(IT));\
            break;\
        }
        switch(itype) {
            CASE_I(0, uint8_t) CASE_I(1, uint16_t) CASE_I(2, uint32_t)
            default: throw std::runtime_error("This should never happen");
        }
        ::close(idfd);
        if(int rc = ::system((std::string("mv ") + idpath + " " + idop).data()); rc != 0) {
            throw std::runtime_error("Failed to mv idpath to final idpath");
        }
    }

    std::fprintf(stderr, "Shape: %zu, %zu;%zu nnz. Max count: %u\n", nlines, njnct, nnz, maxct);

    return 0;
}
