//
// Copyright [] <>
// TODO(srirampc)
//

#ifndef UTILS_HPP
#define UTILS_HPP

#include <cmath>
#include <iostream>
#include <string>
#include <chrono>
#include <ratio>

// Macros for conditional execution
#define PRINT_IF(COND, PRTCD)                                                  \
    {                                                                          \
        if (COND) {                                                            \
            PRTCD;                                                             \
        }                                                                      \
    }


//
template <typename T1, typename T2> struct is_same_type {
    static constexpr bool value = false;
    static constexpr bool yes = false;
    static constexpr bool no = true;
};

template <typename T1> struct is_same_type<T1, T1> {
    static constexpr bool value = true;
    static constexpr bool yes = true;
    static constexpr bool no = false;
};

template <typename eT> inline bool sct_isfinite(eT) { return true; }
template <> inline bool sct_isfinite(float x) { return std::isfinite(x); }
template <> inline bool sct_isfinite(double x) { return std::isfinite(x); }

template <typename eT> inline bool sct_isinf(eT) { return false; }
template <> inline bool sct_isinf(float x) { return std::isinf(x); }
template <> inline bool sct_isinf(double x) { return std::isinf(x); }

class ostream_state {
  private:
    const std::ios::fmtflags orig_flags;
    const std::streamsize orig_precision;
    const std::streamsize orig_width;
    const char orig_fill;

  public:
    explicit inline ostream_state(const std::ostream& o)
        : orig_flags(o.flags()), orig_precision(o.precision()),
          orig_width(o.width()), orig_fill(o.fill()) {}

    inline void restore(std::ostream& o) const {
        o.flags(orig_flags);
        o.precision(orig_precision);
        o.width(orig_width);
        o.fill(orig_fill);
    }
};


template <typename ElemType>
inline void print_elem_zero(const bool modify, std::ostream& o) {

    if (modify) {
        const std::ios::fmtflags save_flags = o.flags();
        const std::streamsize save_precision = o.precision();

        o.unsetf(std::ios::scientific);
        o.setf(std::ios::fixed);
        o.precision(0);

        o << ElemType(0);

        o.flags(save_flags);
        o.precision(save_precision);
    } else {
        o << ElemType(0);
    }
}

template <typename ElemType>
inline void print_elem(const ElemType& x, bool modify, std::ostream& o) {
    if (x == ElemType(0)) {
        return print_elem_zero<ElemType>(modify, o);
    }
    if (std::is_signed<ElemType>::value) {
        if (sct_isfinite(x)) {
            o << x;
        } else {
            o << (sct_isinf(x) ? ((x <= ElemType(0)) ? "-inf" : "inf") : "nan");
        }
    } else {
        o << x;
    }
}

template <typename eT>
inline std::streamsize modify_stream(const eT* data,
                                     const unsigned n_elem, std::ostream& o) {
    o.unsetf(std::ios::showbase);
    o.unsetf(std::ios::uppercase);
    o.unsetf(std::ios::showpos);
    o.fill(' ');

    std::streamsize cell_width;

    bool use_layout_B = false;
    bool use_layout_C = false;
    bool use_layout_D = false;

    for (unsigned i = 0; i < n_elem; ++i) {
        const eT val = data[i];

        if (sct_isfinite(val) == false) {
            continue;
        }

        if (((sizeof(eT) > 4) &&
             (is_same_type<uint64_t, eT>::yes || is_same_type<int64_t, eT>::yes) && 
                 (val >= eT(+10000000000))) ||
            ((sizeof(eT) > 4) &&
             is_same_type<int64_t, eT>::yes && (val <= eT(-10000000000)))) {
            use_layout_D = true;
            break;
        }

        if ((val >= eT(+100)) ||
            ( (std::is_signed<eT>::value) &&
              (val <= eT(-100)) ) ||
            ( (std::is_integral<eT>::value == false) && 
              (val > eT(0)) && (val <= eT(+1e-4)) ) ||
            ( (std::is_integral<eT>::value == false) &&
              (std::is_signed<eT>::value) && 
              (val < eT(0)) && (val >= eT(-1e-4)) )
            //
            // (cond_rel<is_signed<eT>::value>::leq(val, eT(-100))) ||
            // (cond_rel<is_non_integral<eT>::value>::gt(val, eT(0)) &&
            //  cond_rel<is_non_integral<eT>::value>::leq(val, eT(+1e-4))) ||
            // (cond_rel < is_non_integral<eT>::value &&
            //  is_signed<eT>::value > ::lt(val, eT(0)) &&
            //  cond_rel < is_non_integral<eT>::value &&
            //  is_signed<eT>::value > ::geq(val, eT(-1e-4)))
        ) {
            use_layout_C = true;
            break;
        }
        if (
             (val >= eT(+10)) || 
             ( (std::is_signed<eT>::value) && (val <= eT(-10)))
            // (val >= eT(+10)) ||
            // (cond_rel<is_signed<eT>::value>::leq(val, eT(-10)))
        ) {
            use_layout_B = true;
        }
    }

    if (use_layout_D) {
        o.setf(std::ios::scientific);
        o.setf(std::ios::right);
        o.unsetf(std::ios::fixed);
        o.precision(4);
        cell_width = 21;
    } else if (use_layout_C) {
        o.setf(std::ios::scientific);
        o.setf(std::ios::right);
        o.unsetf(std::ios::fixed);
        o.precision(4);
        cell_width = 13;
    } else if (use_layout_B) {
        o.unsetf(std::ios::scientific);
        o.setf(std::ios::right);
        o.setf(std::ios::fixed);
        o.precision(4);
        cell_width = 10;
    } else {
        o.unsetf(std::ios::scientific);
        o.setf(std::ios::right);
        o.setf(std::ios::fixed);
        o.precision(4);
        cell_width = 9;
    }

    return cell_width;
}

template <typename MatType>
inline void eig_print(const MatType& m, const bool modify, std::ostream& o) {

    const ostream_state stream_state(o);

    const std::streamsize cell_width =
        modify ? modify_stream(m.data(), m.size(), o)
               : o.width();

    const unsigned m_n_rows = m.rows();
    const unsigned m_n_cols = m.cols();

    if (m.size() != 0) {
        if (m_n_cols > 0) {
            if (cell_width > 0) {
                for (unsigned row = 0; row < m_n_rows; ++row) {
                    for (unsigned col = 0; col < m_n_cols; ++col) {
                        // the cell width appears to be reset after each element
                        // is printed, hence we need to restore it
                        o.width(cell_width);
                        print_elem(m(row, col), modify, o);
                    }

                    o << '\n';
                }
            } else {
                for (unsigned row = 0; row < m_n_rows; ++row) {
                    for (unsigned col = 0; col < m_n_cols - 1; ++col) {
                        print_elem(m(row, col), modify, o);
                        o << ' ';
                    }

                    print_elem(m(row, m_n_cols - 1), modify, o);
                    o << '\n';
                }
            }
        }
    } else {
        if (modify) {
            o.unsetf(std::ios::showbase);
            o.unsetf(std::ios::uppercase);
            o.unsetf(std::ios::showpos);
            o.setf(std::ios::fixed);
        }

        o << "[matrix size: " << m_n_rows << "x" << m_n_cols << "]\n";
    }

    o.flush();
    stream_state.restore(o);
}

template <typename EigenMatType, typename EigenPrtMatType = EigenMatType>
void eigen_brief_print(const EigenMatType& mx, const std::string header,
                        std::ostream& ostx = std::cout, bool print_size = true) {
    const std::streamsize orig_width = ostx.width();
    if (header.length() != 0) {
        ostx << header << std::endl;
    }
    const ostream_state stream_state(ostx);
    if (print_size) {
        ostx.unsetf(std::ios::showbase);
        ostx.unsetf(std::ios::uppercase);
        ostx.unsetf(std::ios::showpos);
        ostx.setf(std::ios::fixed);
        ostx << "[matrix size: " << mx.rows() << "x" << mx.cols() << "]"
             << std::endl;
    }
    //
    if (mx.size() == 0) {
        ostx.flush();
        stream_state.restore(ostx);
        return;
    }
    //
    if ((mx.rows() <= 5) && (mx.cols() <= 5)) {
        eig_print(mx, true, ostx);
        ostx.flush();
        stream_state.restore(ostx);
        return;
    }
    const bool row_el_flag = (mx.rows() >= 6);
    const bool col_el_flag = (mx.cols() >= 6);
    // Whole matrix
    if (row_el_flag && col_el_flag) {
        EigenPrtMatType px(4, 4);
        px.topLeftCorner(3, 3) = mx.topLeftCorner(3, 3);
        px.topRightCorner(3, 1) = mx.topRightCorner(3, 1);
        px.bottomLeftCorner(1, 3) = mx.bottomLeftCorner(1, 3);
        px.bottomRightCorner(1, 1) = mx.bottomRightCorner(1, 1);
        //
        const std::streamsize cell_width = modify_stream(px.data(), px.size(),
                                                         ostx);
        for (int row = 0; row <= 2; ++row) {
            for(int col=0; col <= 2; ++col) {
                ostx.width(cell_width);
                print_elem(px(row,col), true, ostx);
            }
            ostx.width(6);
            ostx << "...";
            ostx.width(cell_width);
            print_elem(px(row, 3), true, ostx);
            ostx << std::endl;
        }
        for (int col = 0; col <= 2; ++col) {
            ostx.width(cell_width);
            ostx << ":";
        }
        ostx.width(6);
        ostx << "...";
        ostx.width(6);
        ostx << ":" << std::endl;

        const int row = 3; {
            for(int col=0; col <= 2; ++col){
                ostx.width(cell_width);
                print_elem(px(row,col), true, ostx);
            }
            ostx.width(6);
            ostx << "...";
            ostx.width(cell_width);
            print_elem(px(row,3), true, ostx);
            ostx << std::endl;
        }
    }
    //
    if (row_el_flag && !col_el_flag) {
        EigenPrtMatType px(4, mx.cols());
        px.topRows(3) = mx.topRows(3);
        px.bottomRows(1) = mx.bottomRows(1);
        //
        const std::streamsize cell_width = modify_stream(px.data(), px.size(),
                                                         ostx);
        for (int row = 0; row <= 2; ++row) {  // first 3 rows
            for(int col=0; col < mx.cols(); ++col) {
                ostx.width(cell_width);
                print_elem(px(row,col), true, ostx);
            }
            ostx << std::endl;
        }

        for (int col = 0; col < mx.cols(); ++col) {
            ostx.width(cell_width);
            ostx << ":";
        }
        ostx.width(cell_width);
        ostx << std::endl;

        const int row = 3;
        {
            for(int col=0; col < mx.cols(); ++col) {
                ostx.width(cell_width);
                print_elem(px(row,col), true, ostx);
            }
        }
        ostx << std::endl;
    }
    //
    if (!row_el_flag && col_el_flag) {
        EigenPrtMatType px(mx.rows(), 4);
        px.leftCols(3) = mx.leftCols(3);
        px.rightCols(1) = mx.rightCols(1);
        //
        const std::streamsize cell_width = modify_stream(px.data(), px.size(),
                                                         ostx);
        for (int row = 0; row < mx.rows(); ++row) {
            // std::cout << px.row(row).head(3) << "    ..." << px(row, 3)
            //          << std::endl;
            for(int col=0; col <= 2; ++col) {
                ostx.width(cell_width);
                print_elem(px(row,col), true, ostx);
            }
            ostx.width(6);
            ostx << "...";
            ostx.width(cell_width);
            print_elem(px(row,3), true, ostx);
            ostx << std::endl;
        }
    }
    ostx.flush();
    stream_state.restore(ostx);
}

/// macros for block decomposition
#define BLOCK_LOW(i, p, n) ((i * n) / p)
#define BLOCK_HIGH(i, p, n) ((((i + 1) * n) / p) - 1)
#define BLOCK_SIZE(i, p, n) (BLOCK_LOW((i + 1), p, n) - BLOCK_LOW(i, p, n))
#define BLOCK_OWNER(j, p, n) (((p) * ((j) + 1) - 1) / (n))

template <typename SizeType, typename T>
static inline SizeType block_low(const T& rank, const T& nproc,
                                 const SizeType& n) {
    return (rank * n) / nproc;
}

template <typename SizeType, typename T>
static inline SizeType block_high(const T& rank, const T& nproc,
                                  const SizeType& n) {
    return (((rank + 1) * n) / nproc) - 1;
}

template <typename SizeType, typename T>
static inline SizeType block_size(const T& rank, const T& nproc,
                                  const SizeType& n) {
    return block_low<SizeType, T>(rank + 1, nproc, n) -
           block_low<SizeType, T>(rank, nproc, n);
}

template <typename SizeType, typename T>
static inline T block_owner(const SizeType& j, const SizeType& n,
                            const T& nproc) {
    return (((nproc) * ((j) + 1) - 1) / (n));
}

// timer definition
//

template <typename duration> class timer_impl {
  private:
    std::chrono::steady_clock::time_point start;
    typename duration::rep _total_elapsed;

  public:
    const typename duration::rep& total_elapsed = _total_elapsed;

    timer_impl() : start(std::chrono::steady_clock::now()), _total_elapsed(0) {}

    void accumulate() { _total_elapsed += elapsed(); }

    void reset() { start = std::chrono::steady_clock::now(); }

    typename duration::rep elapsed() const {
        std::chrono::steady_clock::time_point stop =
            std::chrono::steady_clock::now();
        typename duration::rep elapsed_time = duration(stop - start).count();
        return elapsed_time;
    }
};

using timer = timer_impl<std::chrono::duration<double, std::milli>>;

#endif  // !UTILS_HPP
