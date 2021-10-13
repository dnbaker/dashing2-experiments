#include <fstream>
#include <iostream>
#include <string>
#include <cassert>
#include <vector>

int main(int argc, char **argv) {
    if(argc != 3) {std::fprintf(stderr, "Expected 3 arguments: executable, inpath, outpath\n"); std::exit(1);}
    std::ios_base::sync_with_stdio(false);
    std::ifstream ifs(argv[1]);
    std::ofstream ofs(argv[2]);
    {
        static constexpr const std::string_view header("#G1\tG2\tK\tsketchsize\tANI\tWJI\tJI\tMash\tDash1\tBD8\tBD4\tBD2\tBD1\tBDN\tSS8\tSS4\tSS2\tSS1\tSSN\tFSS8\tFSS4\tFSS2\tFSS1\tFSSN\tMH8\tMH4\tMH2\tMH1\tMHN\tFMH8\tFMH4\tFMH2\tFMH1\tFMHN\tPMH8Exact\tPMH8-50000000\tPMH8-500000\tPMH4Exact\tPMH4-50000000\tPMH4-500000\tPMH2Exact\tPMH2-50000000\tPMH2-500000\tPMH1Exact\tPMH1-50000000\tPMH1-500000\tPMHNExact\tPMHN-50000000\tPMHN-500000\tBMH8Exact\tBMH8-50000000\tBMH8-500000\tBMH4Exact\tBMH4-50000000\tBMH4-500000\tBMH2Exact\tBMH2-50000000\tBMH2-500000\tBMH1Exact\tBMH1-50000000\tBMH1-500000\tBMHNExact\tBMHN-50000000\tBMHN-500000");
        ofs.write(header.data(), header.size());
        ofs.put('\n');
    }
    std::vector<size_t> offsets;
    for(std::string line; std::getline(ifs, line);ofs.put('\n')) {
        std::string dataline;
        {
            std::ifstream lineifs(line);
            if(!std::getline(lineifs, dataline)) throw 1;
        }
        assert(dataline.size());
        offsets = {size_t(0)};
        for(size_t i = 0; i < dataline.size(); ++i) {
            if(dataline[i] == '\t') offsets.push_back(i);
        }
        ofs.write(dataline.data(), offsets[2] + 1);
        size_t k, sksz;
        auto ptr = &*line.end() - 1;
        for(int ndots = 0; ptr > line.data();--ptr) {
            if(*ptr == '.' && ++ndots == 3) break;
        }
        sksz = std::strtoull(ptr + 1, &ptr, 10);
        k = std::strtoull(ptr + 1, &ptr, 10);
        const std::string dat = std::to_string(k) + "\t" + std::to_string(sksz);
        ofs.write(dat.data(), dat.size());
        ofs.write(&dataline[offsets[2]], dataline.size() - offsets[2]);
    }
}
