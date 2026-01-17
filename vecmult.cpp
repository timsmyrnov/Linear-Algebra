#include <iostream>
#include <vector>
#include <cstdio>
#include <chrono>

struct multiplier {
    using vecd_t = std::vector<double>;
    vecd_t l, r, both;

    multiplier(size_t n) : l(n), r(n), both(2*n) { }

    void init(double shift=1.0) {
        for(size_t i=0; i< l.size(); i++) {
            r[i] = i/1000000.0;
            l[i] = (i+shift)/1000000.0;
            size_t bi = 2*i;
            both[bi] = l[i];
            both[bi+1] = r[i];
        }
    }

    double sep_vec() {
        double res{};
        size_t N = l.size();
        for(size_t i = 0; i< N; ++i) {
            res += l[i]*r[i];
        }
        return res;
    }
    double one_vec() {
        double res{};
        size_t N = both.size();
        for(const double* p = &both[0], *p_end = p+both.size(); p< p_end; p+=2) {
            res += p[0]*p[1];
        }
        return res;
    }
};
auto nanonow() { 
    auto now = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now());
    auto ns_since_unix_epoch = now.time_since_epoch().count();
    return ns_since_unix_epoch;
};

template <typename CB>
uint64_t timeit(CB&& cb) {
    auto t = nanonow();
    cb();
    return nanonow()-t;
}

int main(int argc, char* argv[]) {
    size_t n = (argc> 1? atoi(argv[1]): 1000*1000*10);

    auto t = nanonow();
    std::cerr << "initializing " << n << " elements ...";

    multiplier m(n);
    std::cerr << "done in " << nanonow() -t << std::endl;
    m.init();

    int num_iter = 10;
    uint64_t sep_time_total{}, one_time_total{};
    m.init(0);
    for(int i=0; i< num_iter; i++) {
        double one{};
        double sep{};
        sep_time_total += timeit([&](){ sep = m.sep_vec(); }) ;
        std::cerr << sep_time_total << " " << sep << " vs ";
        m.init(i+2.2);

        one_time_total += timeit([&](){ one = m.one_vec(); }) ;
        std::cerr << one_time_total << " " << one << std::endl;
        m.init(i+1.1);
    }
    std::cerr << "sep to one=" << sep_time_total/(double)one_time_total << std::endl;
}
