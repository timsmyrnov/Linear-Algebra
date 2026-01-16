#pragma once
#include <vector>
#include <iostream>
#include <unordered_map>
#include <optional>
#include <cmath>

namespace matrix {
    template <typename V>
    inline V two_determinant(
        const V& a11, const V& a12,
        const V& a21, const V& a22
    )
    {
        return (a11*a22 - a12*a21);
    }
    template <typename V>
    inline V three_determinant(
        const V& a11, const V& a12, const V& a13, 
        const V& a21, const V& a22, const V& a23,
        const V& a31, const V& a32, const V& a33
    ) 
    {
        return (a11*a22*a33 + a31*a12*a23 + a21*a32*a13 - a31*a22*a13 - a21*a12*a33 - a32*a23*a11);
    }

    namespace square_matrix_exception {
        struct size : public  std::exception {};
        struct numeric : public  std::exception {};
        struct degenerate : public  std::exception {};
    }


    template <typename T>
    class square_matrix {
        using value_t = T;
        using row_t = std::vector<value_t>;
        using data_t = std::vector<row_t>;

        data_t rows;

        size_t size() const { return rows.size(); }

        public:

        // Constructor
        square_matrix(size_t n): rows(n) {
            for(auto& r: rows) { r.resize(n); }
        }
        square_matrix(const data_t& d) : rows(d) {}
        auto& row(size_t i) { return rows[i]; }
        const auto& row(size_t i) const { return rows[i]; }

        auto& get_rows() { return rows; }
        std::ostream& print(std::ostream& os) const {
            os << "[" << std::endl;
            for(const auto& r : rows) {
                os << "[";
                for(const auto& c: r) {
                    os << c << ",";
                }
                os << "]" << std::endl;
            }
            return os << "]";
        }

        template <typename CB>
        void visit_cells(CB&& cb) {
            for(auto& r: rows) {
                for(auto& c: r) {
                    std::forward<CB>(cb)(c);
                }
            }
        }

        template <typename CB>
        void visit_cells(CB&& cb) const {
            for(const auto& r: rows) {
                for(const auto& c: r) {
                    std::forward<CB>(cb)(c);
                }
            }
        }

        template <typename CB>
        void walk(CB&& cb) const {
            for(size_t i =0; i< rows.size(); i++) {
                const auto& r = rows[i];
                for(size_t j =0; j< r.size(); j++) {
                    std::forward<CB>(cb)(r[j], i, j);
                }
            }
        }

        template <typename CB>
        void walk(CB&& cb) {
            for(size_t i = 0; i < rows.size(); i++) {
                auto& r = rows[i];
                for(size_t j = 0; j < r.size(); j++) {
                    std::forward<CB>(cb)(r[j], i, j);
                }
            }
        }

        template <typename CB>
        void diag_walk(CB&& cb) {
            size_t d = size();
            for(size_t i = 0; i < d; i++) {
                std::forward<CB>(cb)(rows[i][i], i);
            }
        }

        template <typename CB>
        void diag_walk(CB&& cb) const {
            size_t d = size();
            for(size_t i = 0; i < d; i++) {
                std::forward<CB>(cb)(rows[i][i], i);
            }
        }

        template <typename CB>
        void col_walk(CB&& cb, int idx) {
            if(idx < 0 || idx >= size()) {
                return;
            }
            for(size_t i = 0; i < size(); i++) {
                std::forward<CB>(cb)(rows[i][idx], i);
            }
        }

        auto& scale(double x) {
            visit_cells([x](auto& v) { v*= x; });
            return *this;
        }

        auto& increment(double x) {
            visit_cells([x](auto& v) { v+= x; });
            return *this;
        }

        auto& replace_col(int idx, const std::vector<double>& col) {
            if (col.size() != size()) {
                return *this;
            }
            col_walk([&](auto& cell, size_t i) { cell = col[i]; }, idx);
            return *this;
        }
        
        struct rowops {
            static void add_mult(row_t& dest, const row_t& src, double factor) {
                for(int i = 0; i < dest.size(); i++) {
                    dest[i] += src[i] * factor;
                }
            }

            static void swap(row_t& x, row_t& y) {
                x.swap(y);
            }

            static void mult(row_t& dest, double factor) {
                dest *= factor;
            }

            static bool cell_is_zero(double x, double epsilon) {
                return std::fabs(x) < std::fabs(epsilon);
            }
        };

        std::pair<bool, size_t> gauss_elimination(double epsilon = EPSILON) {
            size_t N = size(), num_swaps = 0;

            auto fnz = [&](size_t r) -> std::optional<size_t> {
                for(size_t i = r; i < N; i++) {
                    if(!rowops::cell_is_zero(rows[i][r], epsilon))
                        return i;
                }
                return std::nullopt;
            };
            for(int i = 0; i < N; i++) {
                if(auto j = fnz(i); j.has_value()) {
                    if(j.value() > i) {
                        rowops::swap(rows[i], rows[j.value()]);
                        num_swaps++;
                    }
                } else {
                    return std::make_pair(false, num_swaps); // if we have 0's in i-th column for all rows from i-th down - the matrix is degenerate
                }
                row_t& the_row = rows[i];
                for(int j = i + 1; j < N; j++) {
                    row_t& r = rows[j];
                    double the_value = the_row[i];
                    double& ri = r[i];
                    if(!rowops::cell_is_zero(ri, epsilon)) {
                        double factor = - ri / the_value;
                        rowops::add_mult(r, the_row, factor);
                        ri = 0;    
                    }
                }
            }
            return std::make_pair(true, num_swaps);
        }

        value_t is_small() const { return rows.size() < 4; }
        value_t determinant_small() const {
            switch(rows.size()) {
                case 1: return rows[0][0];
                case 2: return two_determinant(
                    rows[0][0], rows[0][1],
                    rows[1][0], rows[1][1]
                );
                case 3: {
                    return three_determinant(
                        rows[0][0], rows[0][1], rows[0][2],
                        rows[1][0], rows[1][1], rows[1][2],
                        rows[2][0], rows[2][1], rows[2][2]
                    );
                }
                default: return {};
            }
        }
        value_t determinant(double epsilon = EPSILON) {
            if(is_small()) return determinant_small();
            else if(auto p = gauss_elimination(epsilon); p.first) {
                double product = 1;
                int mult = 1;
                for(int i = 0; i < p.second; i++) {
                    mult *= -1;
                }

                diag_walk([&product](const auto& c, size_t i){ product *= c; });
                return mult*product;
            }
            return {};
        }
        template <typename CB>
        void zip(CB&& cb, const square_matrix& o) {
            for(size_t i = 0; i < rows.size(); i++) {
                auto& r1 = rows[i];
                auto& r2 = o.rows[i];
                for(size_t j = 0; j < r1.size(); j++) {
                    std::forward<CB>(cb)(r1[j], r2[j], i, j);
                }
            }
        }
        void assert_size(const square_matrix& o) { if(o.size() != size()) {throw square_matrix_exception::size();} }
        auto& operator +=(const square_matrix& o) {
            assert_size(o);
            zip([](auto& a, const auto& b, size_t i, size_t j){ a += b; }, o);
            return *this;
        }
        struct const_col_iterator {
            const square_matrix& m;
            size_t row{};
            const size_t col{}; 
            const_col_iterator () = default;
            const_col_iterator (const const_col_iterator& ) = default;
            const_col_iterator (const square_matrix& x, size_t c, bool is_end=false) : 
                m(x), row(is_end ? m.size() : 0), col(c)
            {}
            const_col_iterator& operator ++() { ++row; return *this; }
            const_col_iterator& operator ++(int) { row++; return *this; }
            const square_matrix::value_t& operator *() const { return m.rows[row][col]; } 
            bool operator !=(const const_col_iterator& o) const { return col != o.col || row != o.row; }
            bool operator ==(const const_col_iterator& o) const { return col == o.col && row == o.row; }
        };

        struct const_row_iterator {
            typename square_matrix::row_t::const_iterator p; 
            const_row_iterator() = default;
            const_row_iterator(const const_row_iterator& ) = default;
            const_row_iterator(const square_matrix& m, size_t i, bool is_end=false) : p(is_end ? m.rows[i].end(): m.rows[i].begin()) {}
            const_row_iterator& operator ++() { ++p; return *this; }
            const_row_iterator& operator ++(int) { p++; return *this; }
            const square_matrix::value_t& operator *() const { return *p; } 
            bool operator !=(const const_row_iterator& o) const {return p != o.p; }
            bool operator ==(const const_row_iterator& o) const {return p == o.p; }
            /* 
            const_row_iterator ri = m.row_begin(i), ri_end = m.row_end(i);
            const_col_iterator ci = m.col_begin(i), ci_end = m.col_end(i);
            for(; ri = row_end(i); ri++, ci++) {
            }
            */
        };
        const_row_iterator row_begin(size_t r) const { return const_row_iterator(*this, r); }
        const_row_iterator row_end(size_t r) const { return const_row_iterator(*this, r, true); }
        std::pair<const_row_iterator, const_row_iterator> row_begin_end(size_t r) const { 
            return std::make_pair(row_begin(r), row_end(r)); 
        }

        const_col_iterator col_begin(size_t r) const { return const_col_iterator(*this, r); }
        const_col_iterator col_end(size_t r) const { return const_col_iterator(*this, r, true); }
        std::pair<const_col_iterator, const_col_iterator> col_begin_end(size_t r) const { 
            return std::make_pair(col_begin(r), col_end(r)); 
        }

        static constexpr double EPSILON = 1e-12;

    };

}
