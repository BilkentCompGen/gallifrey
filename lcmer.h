#ifndef LCMER_H
#define LCMER_H
#include <array>
#include <algorithm>
template<int N>
struct lcmer {
    std::array<uint32_t, N> data;

    lcmer () {}
    lcmer (auto it) {
        std::copy(it, it + N, std::begin(data));
    }

    bool operator==(const lcmer &other) const {
        for (int i = 0; i < N; ++i) {
            if (data[i] != other.data[i]) return false;
        }
        return true;
    }
    bool operator<(const lcmer &other) const {
        for (int i = 0; i < N; ++i) {
            if (data[i] != other.data[i]) {
                if (data[i] < other.data[i]) {
                    return true;
                } else {
                    return false;
                }
            }
        }
        return false;
    }
    void print_name(FILE *out) const {
        fprintf(out, "%0x", data[0]);
        for (int z = 1; z < N; ++z) {
            fprintf(out, "%0x", data[z]);
        }
    }


};
template<int N>
void lcmer_set0(lcmer<N> &l, uint32_t new_core) {
    l.data[N - 1] = new_core;
}
template<int N>
lcmer<N> shift_lcmer(lcmer<N> l) {
    for (int i = 1; i < N; ++i) {
        l.data[i - 1] = l.data[i];
    }
    return l;
}
template<int N>
lcmer<N> update_lcmer(lcmer<N> l, uint32_t new_core) {
    // std::shift_left( std::begin(l.data), std::end(l.data), 1);
    shift_lcmer(l);
    l.data[N - 1] = new_core;
    return l;
}


template <int N>
struct std::hash<lcmer<N>> {
    std::size_t operator()(const lcmer<N> &l) const noexcept {
        return MurmurHash3_32(&l.data, N * sizeof(uint32_t));
    }
};
#endif