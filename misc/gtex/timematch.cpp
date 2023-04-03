#include <string>
#include <cassert>
#include <vector>
#include <chrono>
#include <cstring>

inline bool matchchr(const char *s) {
    // | 32 lower-cases the letter, allowing C and c to match
    return (*s | char(32)) == 'c' && s[1] == 'h' && s[2] == 'r';
}
inline bool matchchr(const std::string &s) {return matchchr(s.data());}

inline bool matchchrman(const char *s) {
    static constexpr uint32_t md = 7497827;
    uint32_t i = 0;
    std::memcpy(&i, s, 3);
    *((uint8_t *)&i) |= 32;
    return i == md;
}
inline bool matchchrman(const std::string &s) {return matchchrman(s.data());}
auto gett() {return std::chrono::high_resolution_clock::now();}

int main() {
    static constexpr size_t n = 1000000;
    std::vector<std::string> rstrings(n);
    std::vector<std::string> cstrings(n);
    for(size_t i = 0; i < n; ++i) {
        rstrings[i] = std::to_string((i * 0x234917123) - 0x14920855l);
        cstrings[i] = std::to_string((i * 0x234917123) - 0x14920855l);
        cstrings[i] = "Chr" + cstrings[i];
#if 0
        if((i & 8) == 0) {
        } else if((i & 16) == 0) {
            cstrings[i] = "chr" + cstrings[i];
            assert(matchchr(cstrings[i]));
        }
#endif
    }
    auto t = gett();
    size_t nmr = 0, cmr = 0;
    for(size_t i = 0; i < n; ++i) {
        nmr += matchchr(rstrings[i]);
    }
    assert(nmr == 0);
    auto t2 = gett();
    for(size_t i = 0; i < n; ++i) {
        cmr += matchchr(cstrings[i]);
    }
    auto t3 = gett();
    assert(cmr > 0);
    auto t1p = std::chrono::duration<double, std::milli>(t2 - t).count();
    auto t2p = std::chrono::duration<double, std::milli>(t3 - t2).count();
    std::fprintf(stderr, "Time for lazy version, never matching: %gms\n", t1p);
    std::fprintf(stderr, "Time for lazy version, sometimes matching: %gms\n", t2p);
    nmr = cmr = 0;
    t = gett();
    for(size_t i = 0; i < n; ++i) {
        nmr += matchchrman(rstrings[i]);
    }
    t2 = gett();
    for(size_t i = 0; i < n; ++i) {
        cmr += matchchrman(cstrings[i]);
    }
    t3 = gett();
    assert(nmr == 0);
    assert(cmr > 0);
    t1p = std::chrono::duration<double, std::milli>(t2 - t).count();
    t2p = std::chrono::duration<double, std::milli>(t3 - t2).count();
    std::fprintf(stderr, "Time for proactive version, never matching: %gms\n", t1p);
    std::fprintf(stderr, "Time for proactive version, sometimes matching: %gms\n", t2p);
}
