#include <mio.hpp>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <iterator>

using std::uint64_t;
using namespace std::string_literals;

void usage() {
    std::fprintf(stderr, "Usage: csrt <flags> data.bin indices.bin indptr.bin remap\n");
    std::fprintf(stderr, "Flags:  \n"
                         "-h: Emit usage\n"
                         "-d <int>: Set (integral) data size in bytes. (Must be 2, 4, or 8). Default: 4\n"
                         "-I <int>: Set (integral) indices size in bytes. (Must be 2, 4, or 8). Default: 4\n"
                         "-P <int>: Set (integral) indptr size in bytes. (Must be 4 or 8). Default: 8\n"
                );
    std::exit(1);
}

static constexpr std::array<int, 3> ISIZE{{2, 4, 8}};

static constexpr bool inisz(int v) {
    return v == 2 || v == 4 || v == 8;
}

size_t getnlines(std::string path) {
    size_t nlines = 0;
    std::ifstream ifs(path);
    for(std::string line;std::getline(ifs, line);) {
        ++nlines;
    }
    return nlines;
}
std::vector<uint64_t, Alloc> &fill_mat(std::vector<uint64_t, Alloc> &r, std::string path, unsigned int ipsize) {
    std::ifstream ifs(path, std::ios::binary);
    if(ipsize == sizeof(uint64_t) {
        r.assign(std::istream_iterator<uint64_t>(ifs), std::istream_iterator<uint64_t>{});
    } else if(ipsize == 2) {
        r.assign(std::istream_iterator<uint16_t>(ifs), std::istream_iterator<uint16_t>{});
    } else if(ipsize == 4) {
        r.assign(std::istream_iterator<uint32_t>(ifs), std::istream_iterator<uint32_t>{});
    } else throw std::runtime_error("fill_mat requires ipsize to be 2, 4, or 8!");
    return r;
}

void perform_tx(const void *dptr, const void *iptr, const std::vector<uint64_t> &indptr, size_t dn, std::vector<std::FILE *> ofps, std::vector<uint64_t> &outindptr, int dsize, int isize) {
    const size_t ninput_rows = indptr.size() - 1;
    for(uint64_t i = 0;i < ninput_rows; ++i) {
        // This corresponds to a feature in the transformed matrix;
        uint64_t ipv, dv;
        for(size_t ij = indptr[i], ej = indptr[i + 1]; ij < ej; ++ij) {
            uint64_t dv;
            if(isize == 2) {
                ipv = ((const uint16_t *)iptr)[ij];
            } else if(isize == 4) {
                ipv = ((const uint32_t *)iptr)[ij];
            } else if(isize == 8) {
                ipv = ((const uint64_t *)iptr)[ij];
            } else {throw 5;}
            if(dsize == 2) {
                dv = ((const uint16_t *)dptr)[ij];
            } else if(isize == 4) {
                dv = ((const uint32_t *)dptr)[ij];
            } else if(isize == 8) {
                dv = ((const uint64_t *)dptr)[ij];
            } else {throw 5;}
            auto fp = ofps.at(ipv);
            if(std::fwrite(&i, sizeof(i), 1, fp) != 1) throw 4;
        }
    }
    std::vector<uint64_t> outindptr(ofps.size() + 1);
    // Since we want the first number to be zero, we allocate one extra integer, and do the partial sum with the 0 included.
    std::transform(ofps.data(), &ofps[ofps.size()], outindptr.begin() + 1, [](std::FILE *fp) {
        const int fd = ::fileno(fp);
        struct stat st;
        if(::fstat(fd, &st)) {std::perror("Failed to stat"); std::exit(1);}
        return st.st_size / sizeof(uint64_t);
    });
    std::partial_sum(outindptr.begin(), outindptr.end(), outindptr.begin());
}


#define manual_assert(x) do {if(__builtin_expect(!(x), 0)) {std::cerr << #x << '\n'; throw 1;} } while(0)



int main(int argc, char **argv) {
    int dsize = 4, isize = 4, ipsize = 8;
    for(int c;(c = getopt(argc, argv, "h?")) >= 0;) { switch(c) {
        default: case 'h': case '?': usage(); break;
        case 'd': dsize = std::atoi(optarg); if(!inisz(dsize)) goto FAIL; break;
        case 'I': isize = std::atoi(optarg); if(!inisz(isize)) goto FAIL; break;
        case 'P': ipsize = std::atoi(optarg); if(!inisz(ipsize)) goto FAIL; break;
    }}
    if(0) {
        FAIL:
        std::fprintf(stderr, "Sizes must be 2, 4 or 8\n");
        usage();
    }
    auto diff = argc - optind;
    if(diff != 4) {
        std::fprintf(stderr, "Expected 4 positional arguments\n");
        usage();
    }
    uint64_t rv = 0;
    for(auto p = argv; *p; ++p) rv ^= std::hash<std::string>()(*p);

    std::string dpath = argv[optind], ipath = argv[optind + 1], ippath = argv[optind + 2], rmpath = argv[optind + 3];
    std::vector<uint64_t> input_indptrs;
    fill_mat(input_indptrs, ippath, ipsize);
    if(input_indptrs.empty()) throw 1;
    const size_t input_nrows = input_indptrs.size() - 1;
    auto tmpdir = "tmp."s + std::to_string(rv);
    ::system(("mkdir "s + tmpdir).data());
    const size_t input_ncols = getnlines(rmpath);
    std::vector<std::FILE *> ofpds(input_ncols);
    std::vector<std::FILE *> ofpis(input_ncols);
    for(size_t i = 0; i < input_ncols; ++i) {
        std::string path = tmpdir + "/" + std::to_string(i) + + ".data.bin";
        if((ofpds[i] = std::fopen(path.data(), "wb")) == nullptr)
            throw std::runtime_error("Failed to open "s + path + " for writing");
        path = tmpdir + "/" + std::to_string(i) + + ".id.bin";
        if((ofpis[i] = std::fopen(path.data(), "wb")) == nullptr)
            throw std::runtime_error("Failed to open "s + path + " for writing");
    }
    mio::mmap_sink data_mmap(dpath);
    mio::mmap_sink indices_mmap(ipath);
    const void *dptr = static_cast<const void *>(data_mmap.data());
    const void *iptr = static_cast<const void *>(indices_mmap.data());
    std::vector<uint64_t> output_indptrs;
    perform_tx(dptr, iptr, input_indptrs, dn, ofpds, ofpis, output_indptrs);
    for(auto fp: ofpds) std::fclose(fp);
    for(auto fp: ofpis) std::fclose(fp);
    ::system(("rm -r "s + tmpdir).data());
}
