#pragma once
#include <vector>
#include <iostream>

namespace linear_algebra {

    template <typename T>
    class matrix {
        using value_t = T;
        using row_t = std::vector<value_t>;
        using data_t = std::vector<row_t>;

        data_t rows;

        size_t height() const { return rows.size(); }
        size_t width() const { return (rows.empty() ? 0: rows[0].size()); }
        size_t size() const { return std::min(height(), width()); }
        bool is_square() const { return height() == width(); }

        public:

        // Constructor
        matrix(size_t n, size_t m): rows(n) {
            for(auto& r: rows) { r.resize(m); }
        }
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
            if (idx < 0 || idx >= width()) {
                return;
            }
            size_t m_height = height();
            for(size_t i = 0; i < m_height; i++) {
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

        auto& insert_col(int idx, std::vector<double> col) {
            size_t m_height = height();
            size_t c_size = col.size();
            if (c_size != m_height) {
                return *this;
            }
            col_walk([&](auto& cell, size_t i) { cell = col[i]; }, idx);
            return *this;
        }
    };
}
