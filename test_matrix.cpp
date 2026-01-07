#include "include/matrix.hpp"
#include <vector>
#include <unordered_map>
#include <optional>

using sqmatrix_double_t = matrix::square_matrix<double>;

int basic_test(int argc, char* argv[]) {
    matrix::square_matrix<double> mtx(5);
    std::vector<double> column = {0, 0, 0, 0, 0};

    { // Direct access
    double i = 0;
    for(auto& r : mtx.get_rows()) {
        for(auto& c: r) {
            c = i++ / 10;
        }
    }
    mtx.print(std::cout) << std::endl;
    }

    { // Visitor access
        double x = 0.1;
        size_t count = 0;
        mtx.visit_cells([&x, &count](auto& c){
            c += x;
            x += 0.01;
            count++;
        });
        mtx.print(std::cout) << std::endl;
        std::cout << "COUNT=" << count << std::endl;
    }

    { // Fancy visitor
        struct fancy_cb {
            size_t c{};
            double increment{};
            fancy_cb(double i) : increment(i){}

            void operator() (double& x) {
                c++;
                x += increment;
            }
        };

        fancy_cb cb(0.3);

        mtx.visit_cells(cb);
        mtx.visit_cells([](const auto& c){ std::cout << c << ","; });
        std::cout <<  std::endl << " COUNT:" << cb.c << std::endl;
        mtx.visit_cells(cb);
        std::cout <<  " COUNT:" << cb.c << std::endl;
    }

    {
        mtx.visit_cells([](auto& c){ c = 0; });
        mtx.diag_walk([](auto& c, size_t i){ c = 1; });
        mtx.print(std::cout) << std::endl;
        mtx.scale(10).increment(5).print(std::cout) << std::endl;
        mtx.replace_col(3, column);
        mtx.print(std::cout) << std::endl;
    }
    return 0;
}

namespace {
    void test_det(sqmatrix_double_t&& m) {
        m.print(std::cout << "BEFORE: ");
        bool elim_status = m.gauss_elimination();
        m.print(std::cout << "AFTER: ") << (elim_status? "SUCCESS": "FAIL") <<std::endl;
    }
}
int det_test(int argc, char* argv[]) {
    test_det(sqmatrix_double_t{{{7,1}, {6,2}}}); // success
    test_det(sqmatrix_double_t{{{4,6}, {2,3}}}); // fail
    test_det(sqmatrix_double_t{{{7,1,3}, {14,2,6}, {2,3,5}}}); // success
    test_det(sqmatrix_double_t{{{7,1,3}, {1,2,6}, {2,3,5}}}); // fail
    test_det(sqmatrix_double_t{{{1,1,1}, {2,3,3}, {2,2,3}}}); // success
    test_det(sqmatrix_double_t{{{1,0,0}, {0,1,0}, {0,0,1}}}); // success
    return 0;
}

int main(int argc, char* argv[]) {
    // basic_test(argc, argv);
    det_test(argc, argv);
}