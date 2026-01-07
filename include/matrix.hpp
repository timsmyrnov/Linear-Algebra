#pragma once
#include <vector>
#include <iostream>
#include <unordered_map>
#include <optional>
#include <cmath>

namespace matrix {

    struct submatrix {
        using bitvec_t = std::vector<bool>;
        bitvec_t rows, cols; // true means included 
        void exclude(size_t i, size_t j) {
            rows[i] = cols[j] = false;
        }
        submatrix(size_t n) : rows(n), cols(n) {}
    };

    template <typename T>
    struct submatrix_hash {
        using value_t = T;
        using colmap_t = std::unordered_map<submatrix::bitvec_t, value_t>;
        using rowmap_t = std::unordered_map<submatrix::bitvec_t, colmap_t>;
        rowmap_t map;

        std::optional<value_t> get_value(const submatrix::bitvec_t& rows, const submatrix::bitvec_t& cols) const {
            if(auto i = map.find(rows); i != map.end()) {
                if(auto j = i->second.find(cols); j != i->second.end())
                    return j->second;
            }
            return std::nullopt;
        } 
        void set_value(const submatrix::bitvec_t& rows, const submatrix::bitvec_t& cols, const value_t& v) {
            typename rowmap_t::iterator rmi = map.find(rows);
            
            if(rmi == map.end())
                rmi = map.insert(std::make_pair(rows, colmap_t{})).first;

            colmap_t& colmap = *(rmi->second);
            typename colmap_t::iterator cmi = colmap.find(cols);
            if(cmi == colmap.end())
                cmi = colmap.insert(std::make_pair(cols, value_t{})).first;
                
            cmi->second = v;
        }
    };

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

        bool gauss_elimination(double epsilon = 1e-12) {
            size_t N = size();
            auto fnz = [&](size_t r) -> std::optional<size_t> {
                for(size_t i = r; i < N; i++) {
                    if(!rowops::cell_is_zero(rows[i][r], epsilon))
                        return i;
                }
                return std::nullopt;
            };
            for(int i = 0; i < N; i++) {
                if(auto j = fnz(i); j.has_value()) {
                    if(j.value() > i)
                        rowops::swap(rows[i], rows[j.value()]);
                } else {
                    return false; // if we have 0's in i-th column for all rows from i-th down - the matrix is degenerate
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
            return true;
        }
    };
}
