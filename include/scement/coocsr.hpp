//
// Copyright [2024] <srirampc>
//

#ifndef COOCSR_HPP
#define COOCSR_HPP

#include <cassert>
#include <utility>
#include <vector>
// OpenMP
#include "scement/utils.hpp"
#include <omp.h>

template <typename fmT, typename imT, typename eT = float,
          typename iT = uint32_t>
class COOCSR {
  public:
    using FMatType = fmT;
    using IMatType = imT;
    using ElemType = eT;
    using IndexType = iT;

  private:  // Size
    IndexType _n_row;
    IndexType _n_cols;
    IndexType _n_elem;
    FMatType _data;
    IMatType _coo;
    IMatType _csrptr;
    // IMatType _csrindex;
    bool load_flag;

  public:  // Public properties
    const IndexType& n_rows = _n_row;
    const IndexType& n_cols = _n_cols;
    const IndexType& n_elem = _n_elem;
    const FMatType& data = _data;
    const IMatType& coo = _coo;
    const IMatType& csrptr = _csrptr;

  public:
    COOCSR() {
        _n_row = _n_cols = _n_elem = 0;
        load_flag = false;
    }

    // COOCSR(IndexType nr, IndexType nc, IndexType nel, const FMatType& dat,
    //        const IMatType& icoo, const IMatType& irptr)
    //     : _n_row(nr), _n_cols(nc), _n_elem(nel), _data(dat), _coo(icoo),
    //       _csrptr(irptr), load_flag(true) {}

    COOCSR(IndexType nr, IndexType nc, IndexType nel, FMatType&& dat,
           IMatType&& icoo, IMatType&& irptr)
        : _n_row(nr), _n_cols(nc), _n_elem(nel), _data(std::move(dat)),
          _coo(std::move(icoo)), _csrptr(std::move(irptr)), n_rows(_n_row),
          n_cols(_n_cols), n_elem(_n_elem), data(_data), coo(_coo),
          csrptr(_csrptr), load_flag(true) {
        // Update References
        //
        // _n_row = nr;
        // _n_cols = nc;
        // _n_elem = nel;
        // data = _data;
        // coo = _coo;
        // csrptr = _csrptr;
    }

    // COOCSR(const COOCSR& rcx)
    //     : _n_row(rcx._n_row), _n_cols(rcx._n_cols), _n_elem(rcx._n_elem),
    //       _data(rcx._data), _coo(rcx._coo), _csrptr(rcx._csrptr),
    //       load_flag(rcx.load_flag) {}

    COOCSR(COOCSR&& rcx)
        : _n_row(rcx._n_row), _n_cols(rcx._n_cols), _n_elem(rcx._n_elem),
          _data(std::move(rcx._data)), _coo(std::move(rcx._coo)),
          _csrptr(std::move(rcx._csrptr)), n_rows(_n_row), n_cols(_n_cols),
          n_elem(_n_elem), data(_data), coo(_coo), csrptr(_csrptr),
          load_flag(true) {
        // _n_row = rcx.n_rows;
        // _n_cols = rcx.n_cols;
        // _n_elem = rcx.n_elem;
        // n_rows = _n_row;
        // n_cols = _n_cols;
        // n_elem = _n_elem;
        // data = _data;
        // coo = _coo;
        // csrptr = _csrptr;
    }

    COOCSR& operator=(COOCSR&& rcx) {
        _n_row = rcx._n_row;
        _n_cols = rcx._n_cols;
        _n_elem = rcx._n_elem;
        _data = std::move(rcx._data);
        _coo = std::move(rcx._coo);
        _csrptr = std::move(rcx._csrptr);
        return *this;
    }

    std::size_t size() const { return _n_elem; }

    void set_data_at(std::size_t index, ElemType value) {
        assert(index < n_elem);
        _data[index] = value;
    }

    template <typename RowSelType, typename ColSelType>
    void set(const RowSelType& row_select, const ColSelType& col_select,
             ElemType value) {
        assert(row_select.size() == _n_row);
        assert(col_select.size() == _n_cols);
#pragma omp parallel for
        for (std::size_t i = 0; i < _n_elem; i++) {
            std::size_t row_idx = _coo(i, 0);
            std::size_t col_idx = _coo(i, 1);
            if (col_select(col_idx) && row_select(row_idx)) {
                _data(i) = value;
            }
        }
    }

};

template <typename COOCSRType,
          typename FMatType = typename COOCSRType::FMatType>
class COOCSRRef {
  public:
    using ElemType = typename COOCSRType::ElemType;
    using IndexType = typename COOCSRType::IndexType;
    using IMatType = typename COOCSRType::IMatType;

  private:
    FMatType _data;

  public:
    const std::size_t n_rows;
    const std::size_t n_cols;
    const std::size_t n_elem;
    const FMatType& data;
    const IMatType& coo;
    const IMatType& csrptr;

    explicit COOCSRRef(const COOCSRType& src, const FMatType& dx)
        : n_rows(src.n_rows), n_cols(src.n_cols), n_elem(src.n_elem), _data(dx),
          data(_data), coo(src.coo), csrptr(src.csrptr) {}

    void set_data_at(std::size_t index, ElemType value) {
        assert(index < n_elem);
        _data(index) = value;
    }

    template <typename RowSelType, typename ColSelType>
    void set(const RowSelType& row_select, const ColSelType& col_select,
             ElemType value) {
        assert(row_select.size() == n_rows);
        assert(col_select.size() == n_cols);
#pragma omp parallel for
        for (std::size_t i = 0; i < n_elem; i++) {
            std::size_t row_idx = coo(i, 0);
            std::size_t col_idx = coo(i, 1);
            if (col_select(col_idx) && row_select(row_idx)) {
                _data(i) = value;
            }
        }
    }
};

template <typename imT, typename cvT, typename ivT, typename iT = uint32_t>
class SparseCfg {
  public:
    using IMatType = imT;
    using CVecType = cvT;
    using LVecType = ivT;
    using IndexType = iT;

  private:
    IMatType _unique;
    CVecType _counts;
    LVecType _lookup;
    IndexType _n_config;
    bool load_flag;

  public:
    const IMatType& unique = _unique;
    const CVecType& counts = _counts;
    const LVecType& lookup = _lookup;
    const unsigned& n_config = _n_config;

  public:
    SparseCfg() { load_flag = false; }
    // SparseCfg(const IMatType& unq, const CVecType& ctx, const LVecType& lkp,
    //           const IndexType ncfg)
    //     : _unique(unq), _counts(ctx), _lookup(lkp), _n_config(ncfg) {
    //     load_flag = true;
    // }

    SparseCfg(IMatType&& unq, CVecType&& ctx, LVecType&& lkp, IndexType ncfg)
        : _unique(std::move(unq)), _counts(std::move(ctx)),
          _lookup(std::move(lkp)), _n_config(ncfg), unique(_unique),
          lookup(_lookup), counts(_counts), n_config(_n_config),
          load_flag(true) {
        // unique = _unique;
        // counts = _counts;
        // lookup = _lookup;
        // _n_config = ncfg;
        // load_flag = true;
    }

    // SparseCfg(const SparseCfg& rcx)
    //     : _n_config(rcx._n_config), _unique(rcx._unique),
    //     _lookup(rcx._lookup),
    //       _counts(rcx._counts) {}

    SparseCfg(SparseCfg&& rcx)
        : _unique(std::move(rcx._unique)), _lookup(std::move(rcx._lookup)),
          _counts(std::move(rcx._counts)), _n_config(rcx._n_config),
          unique(_unique), lookup(_lookup), counts(_counts),
          n_config(_n_config), load_flag(true) {
        // unique = _unique;
        // counts = _counts;
        // lookup = _lookup;
        // _n_config = rcx.n_config;
        // load_flag = true;
    }

    SparseCfg& operator=(SparseCfg&& rcx) {
        _n_config = rcx._n_config;
        _unique = std::move(rcx._unique);
        _lookup = std::move(rcx._lookup);
        _counts = std::move(rcx._counts);
        // update refs
        load_flag = rcx.load_flag;
        return *this;
    }

    friend class HDF5LoaderA;
    friend class HDF5LoaderE;
};

template <typename SparseCfgType, typename SelectType> struct SparseCfgSlice {
  public:
    using IMatType = typename SparseCfgType::IMatType;
    using CVecType = typename SparseCfgType::CVecType;
    using LVecType = typename SparseCfgType::LVecType;
    using IndexType = typename SparseCfgType::IndexType;

  public:
    const IMatType& unique;
    const CVecType& counts;
    const LVecType& lookup;
    const SelectType& select;
    const unsigned n_config;

  public:
    //
    explicit SparseCfgSlice(const IMatType& unq, const CVecType& cts,
                            const LVecType& lkp, const unsigned& ncfg,
                            const SelectType& sel)
        : unique(unq), counts(cts), lookup(lkp), select(sel), n_config(ncfg) {}
    //
    explicit SparseCfgSlice(const SparseCfgType& cfg, const SelectType& sel)
        : unique(cfg.unique), counts(cfg.counts), lookup(cfg.lookup),
          select(sel), n_config(cfg.n_config) {}
};

//
// Element-wise coocsr dense multiplication
//
template <typename InSpMatType, typename DenseVecType, typename OutSpMatType>
void elt_mult_coocsr_dense(const InSpMatType& spmat, const DenseVecType& fdense,
                           OutSpMatType* result) {
    assert(spmat.n_rows == fdense.size());
    assert(spmat.n_elem == result->n_elem);
    using ElemType = typename InSpMatType::ElemType;
    //
#pragma omp parallel for
    for (std::size_t k = 0; k < spmat.n_elem; k++) {
        auto row_idx = spmat.coo(k, 0);
        ElemType sval = spmat.data(k),
                 dval = static_cast<ElemType>(fdense(row_idx));
        result->set_data_at(k, sval * dval);
    }
}

template <typename OutDenseType, typename OutElemType, typename InDenseType,
          typename SpMatType>
OutDenseType mult_coocsr_dense(const SpMatType& spmat,
                               const InDenseType& fdense, std::size_t n_outrows,
                               std::size_t n_outcols) {
    assert(spmat.n_rows == n_outrows);
    using IndexType = typename SpMatType::IndexType;
    //
    // Distribute Rows
    //
    int nthreads = omp_get_max_threads();
    std::vector<IndexType> thread_starts(nthreads + 1, 0);
    for (auto i = 1; i < nthreads; i++) {
        auto bstart = block_low(i, nthreads, spmat.n_elem);
        thread_starts[i] = spmat.coo(bstart, 0) + 1;
    }
    thread_starts[nthreads] = spmat.n_rows;
    // thread_starts.print("TSTART:");
    //
    OutDenseType result(n_outrows, n_outcols);
#pragma omp parallel for
    for (IndexType idx = 0; idx < result.size(); idx++) {
        result(idx) = 0;
    }

#pragma omp parallel for
    for (int t = 0; t < nthreads; t++) {
        for (IndexType srow = thread_starts[t]; srow < thread_starts[t + 1];
             srow++) {
            IndexType begin_idx = spmat.csrptr(srow, 0);
            IndexType end_idx = spmat.csrptr(srow, 1);
            for (IndexType sidx = begin_idx; sidx < end_idx; sidx++) {
                IndexType scol = spmat.coo(sidx, 1);
                OutElemType sval = static_cast<OutElemType>(spmat.data(sidx));
                for (IndexType dcol = 0; dcol < n_outcols; dcol++) {
                    OutElemType dval =
                        static_cast<OutElemType>(fdense(scol, dcol));
                    OutElemType uval = result(srow, dcol) + (sval * dval);
                    result(srow, dcol) = uval;
                }
            }
        }
    }
    //
    return result;
}

template <typename DenseVecType, typename SpMatType, typename DenseMatType,
          typename CfgType>
DenseVecType coocsr_row_sum(const SpMatType& xdata, const DenseMatType& ydata,
                            const CfgType& config) {
    // assert(xdata.n_rows == ydata.n_rows);          // TODO(me)
    // assert(config.n_config == ydata.n_cols);       // TODO(me)
    // assert(config.counts.n_rows == ydata.n_cols);  // TODO(me)
    assert(xdata.n_cols == config.select.size());
    using ElemType = typename SpMatType::ElemType;
    using IndexType = typename SpMatType::IndexType;
    //
    // Distribute Rows
    //
    int nthreads = omp_get_max_threads();
    std::vector<IndexType> thread_starts(nthreads + 1, 0);
    for (int i = 1; i < nthreads; i++) {
        IndexType bstart = block_low(i, nthreads, xdata.n_elem);
        thread_starts[i] = xdata.coo(bstart, 0) + 1;
    }
    thread_starts[nthreads] = xdata.n_rows;
    // thread_starts.print("TSTART:");
    // rsum
    DenseMatType rsum(xdata.n_rows, config.n_config);
#pragma omp parallel for
    for (IndexType i = 0; i < rsum.size(); i++) {
        rsum(i) = ElemType(0);
    }
    //
#pragma omp parallel for
    for (int t = 0; t < nthreads; t++) {
        for (IndexType row_id = thread_starts[t]; row_id < thread_starts[t + 1];
             row_id++) {
            IndexType begin_idx = xdata.csrptr(row_id, 0);
            IndexType end_idx = xdata.csrptr(row_id, 1);
            for (IndexType sidx = begin_idx; sidx < end_idx; sidx++) {
                IndexType col_id = xdata.coo(sidx, 1);
                IndexType cfg_id = config.lookup(col_id);
                if (config.select(col_id) == 0)
                    continue;
                rsum(row_id, cfg_id) += xdata.data(sidx);
            }
        }
    }

    // ysum
#pragma omp parallel for
    for (IndexType col = 0; col < config.n_config; col++) {
#pragma omp parallel for
        for (IndexType row = 0; row < xdata.n_rows; row++) {
            // rsum.at(row, col) -= ysum.at(row, col);
            rsum(row, col) -=
                (ydata(row, col) * static_cast<ElemType>(config.counts(col)));
        }
    }
    // rsum -= ysum;
    DenseVecType rsx(xdata.n_rows);
    // arma::sum(rsum, 1);
#pragma omp parallel for
    for (IndexType row = 0; row < xdata.n_rows; row++) {
        rsx(row) = 0;
        for (IndexType col = 0; col < config.n_config; col++) {
            rsx(row) += rsum(row, col);
        }
    }
    // rsx.head_rows(3).print("rsx");
    return rsx;
}

template <typename DenseVecType, typename SpMatType, typename DenseMatType,
          typename CfgType>
DenseVecType coocsr_row_mean(const SpMatType& xdata, const DenseMatType& ydata,
                             const CfgType& config) {
    using ElemType = typename SpMatType::ElemType;
    using IndexType = typename SpMatType::IndexType;
    ElemType rtotal = static_cast<ElemType>(0);
    //
#pragma omp parallel for reduction(+ : rtotal)
    for (IndexType col = 0; col < config.select.size(); col++) {
        ElemType indf = static_cast<ElemType>(config.select[col]);
        rtotal += indf;
    }

    DenseVecType ravg =
        coocsr_row_sum<DenseVecType, SpMatType, DenseMatType, CfgType>(
            xdata, ydata, config);
    ravg /= rtotal;
    return ravg;
}

template <typename DenseVecType, typename SpMatType, typename DenseMatType,
          typename CfgType>
DenseVecType row_sum_of_squares(const SpMatType& xdata,
                                const DenseMatType& ydata,
                                const CfgType& config) {
    // assert(xdata.n_rows == ydata.n_rows);
    // assert(config.n_config == ydata.n_cols);
    assert(xdata.n_cols == config.lookup.size());
    assert(xdata.n_cols == config.select.size());
    using ElemType = typename SpMatType::ElemType;
    using IndexType = typename SpMatType::IndexType;
    //
    // Distribute Rows
    //
    int nthreads = omp_get_max_threads();
    std::vector<IndexType> thread_starts(nthreads + 1, 0);
    for (auto i = 1; i < nthreads; i++) {
        auto bstart = block_low(i, nthreads, xdata.n_elem);
        thread_starts[i] = xdata.coo(bstart, 0) + 1;
    }
    thread_starts[nthreads] = xdata.n_rows;
    //
    ElemType nselect = 0.0;
#pragma omp parallel for reduction(+ : nselect)
    for (auto i = 0; i < config.select.size(); i++) {
        nselect += static_cast<ElemType>(config.select[i]);
    }
    //
    DenseVecType result(xdata.n_rows);
#pragma omp parallel for
    for (IndexType i = 0; i < result.size(); i++) {
        result(i) = ElemType(0);
    }
    //
#pragma omp parallel for
    for (int t = 0; t < nthreads; t++) {
        for (IndexType trx = thread_starts[t]; trx < thread_starts[t + 1];
             trx++) {
            IndexType begin_idx = xdata.csrptr(trx, 0);
            IndexType end_idx = xdata.csrptr(trx, 1);
            for (IndexType rx = begin_idx; rx < end_idx; rx++) {
                IndexType col = xdata.coo(rx, 1);
                if (!config.select(col))
                    continue;
                IndexType row = xdata.coo(rx, 0);
                IndexType cfg = config.lookup(col);
                ElemType val = xdata.data(rx) - (2 * ydata(row, cfg));
                result(row) += xdata.data(rx) * val;
            }
        }
    }
    //
    // ysq
    DenseMatType cfbsq(ydata);
#pragma omp parallel for
    for (IndexType cx = 0; cx < config.n_config; cx++) {
#pragma omp parallel for
        for (IndexType row = 0; row < xdata.n_rows; row++) {
            cfbsq(row, cx) = static_cast<float>(config.counts(cx)) *
                             cfbsq(row, cx) * cfbsq(row, cx);
        }
    }
    //
    DenseVecType ydsq(xdata.n_rows);
#pragma omp parallel for
    for (IndexType row = 0; row < ydsq.size(); row++) {
        ydsq(row) = 0;
        for (IndexType col = 0; col < config.n_config; col++) {
            ydsq(row) += cfbsq(row, col);
        }
    }
    return result + ydsq;
}

#endif  // !COOCSR_HPP
