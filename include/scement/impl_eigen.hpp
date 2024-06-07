//
// Copyright [2024]
//
#ifndef SCEMENT_EIGEN_IMPL_HPP
#define SCEMENT_EIGEN_IMPL_HPP

#include "scement/coocsr.hpp"
#include "scement/data_if.hpp"
#include "scement/impl_base.hpp"
#include "scement/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <string>

class ScementEigenInterface {
  public:
    //
    using ElemType = float;
    using IntType = int32_t;
    using UIntType = uint32_t;
    using IndexType = uint32_t;
    using SelectType = uint8_t;
    //
    using FMatType = Eigen::MatrixXf;
    using FVecType = Eigen::VectorXf;
    using FRowType = Eigen::Matrix<float, 1, Eigen::Dynamic>;
    using IMatType = Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic>;
    using IVecType = Eigen::Matrix<int32_t, Eigen::Dynamic, 1>;
    using UMatType = Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic>;
    //
    using MapFMatType = Eigen::Map<FMatType>;
    using MapIMatType = Eigen::Map<IMatType>;
    using MapFVecType = Eigen::Map<FVecType>;
    using MapIVecType = Eigen::Map<IVecType>;
    //
    static inline FMatType transpose(const FMatType& infm) {
        return infm.transpose();
    }
    static inline IMatType transpose(const IMatType& infm) {
        return infm.transpose();
    }

    static inline ElemType pinv_tolerance() { return 1e-8; }
    static inline FMatType imat2fmat(const IMatType& infm) {
        return infm.cast<ElemType>();
    }
    static inline FMatType pseudo_inverse(const FMatType& infm) {
        return infm.completeOrthogonalDecomposition().pseudoInverse();
    }
    static inline FMatType head_columns(const FMatType& infm,
                                        std::size_t ncols) {
        return infm.leftCols(ncols);
    }
    static inline FMatType diag(const FRowType& rx) { return rx.asDiagonal(); }
    static inline FRowType colsum(const FMatType& infm) {
        return infm.colwise().sum();
    }
    static inline FVecType colmean(const FMatType& infm) {
        return infm.colwise().mean();
    }
    static inline FVecType colvar(const FMatType& infm) {
        FRowType rx = colmean(infm);
        FMatType rw = (infm.rowwise() - rx).array().square();
        return rw.colwise().sum() / (rw.rows() - 1);
    }
    static inline IntType sum(const IMatType& infm) { return infm.sum(); }
    static inline ElemType mean(const FMatType& infm) { return infm.mean(); }
    static inline ElemType var(const FMatType& infm) {
        return (infm.array() - infm.mean()).square().sum() / (infm.size() - 1);
        // return arma::var(infm.as_col());
    }

    // Key functions
    static inline FMatType gram_matrix(const FMatType& infm) {
        return infm.transpose() * infm;
    }
    static inline FMatType bhat(const FMatType& des_adt,
                                const FMatType& des_gram) {
        // FMatType des_pv = pseudo_inverse(des_gram);  // NCVR x NCVR
        // _B_hat = des_adt * transpose(des_pv);        // NFEATURESxNCVR
        return des_adt * des_gram.completeOrthogonalDecomposition()
                             .pseudoInverse()
                             .transpose();
    }
    static inline FMatType cfg_bhat(const FMatType& bhx,
                                    const IMatType& cfg_uq) {
        // _B_hat *
        //  ScementData::imat2fmat(ScementData::transpose(_config.unique));
        return bhx * cfg_uq.transpose().cast<ElemType>();
    }
    static inline FMatType grand_mean(const FMatType& bhx,
                                      const FMatType& ctx) {
        // _grand_mean =
        //     ScementData::head_columns(_B_hat, _nbatches) *
        //   (ScementData::transpose(_batch_counts) / _narray);  // NFEATURESx1
        // Num of total entries
        // ElemType narray = 0.0;
        // for (IndexType i = 0; i < bt_counts.size(); i++) {
        //     narray += bt_counts(i);
        // }
        return bhx.leftCols(ctx.size()) * (ctx / ctx.sum());
    }
    static inline FMatType ghat_mprime(const FMatType& mp,
                                       const IMatType& ctx) {
        return mp * ctx.cast<ElemType>();
        //_m_prime * SDT::imat2fmat(_ebt_counts);  // NFEATURESxNBATCH
    }

    static inline FMatType ghat_base(const FMatType& bdes) {
        return bdes * FMatType(bdes.colwise().sum().asDiagonal())
                          .completeOrthogonalDecomposition()
                          .pseudoInverse();
    }
};

template <typename EGT, typename EIMT> class ScementEigenPrintHelper {

  public:
    template <typename EigenMatType, typename EigenPrtMatType = EigenMatType>
    static void brief_print(const EigenMatType& mx, const std::string header,
                            std::ostream& ox) {
        eigen_brief_print<EigenMatType, EigenPrtMatType>(mx, header, ox);
    }

    static void print(const typename EGT::IMatType& xmat,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        brief_print(xmat, prt_prefix, ox);
    }

    static void print(const typename EGT::FMatType& xmat,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        brief_print(xmat, prt_prefix, ox);
    }

    static void print(const typename EGT::MapIMatType& xmat,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        brief_print(xmat, prt_prefix, ox);
    }

    static void print(const typename EGT::MapFMatType& xmat,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        brief_print(xmat, prt_prefix, ox);
    }

    static void print(const typename EGT::FVecType& xmat,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        brief_print(xmat, prt_prefix, ox);
    }

    static void print(const typename EGT::SpMatType& spmat,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        // TODO(me) outsource the printer
        ox << prt_prefix << std::endl
           << "->SIZE : [" << spmat.n_rows << " x " << spmat.n_cols << " : "
           << spmat.n_elem << "]" << std::endl;
        brief_print<typename EGT::SpMatType::IMatType, typename EGT::IMatType>(
            spmat.coo, "->COO", ox);
        brief_print<typename EGT::SpMatType::IMatType, typename EGT::IMatType>(
            spmat.csrptr, "->CSRPTR", ox);

        brief_print(spmat.data, "->DATA", ox);
    }

    static void print(const typename EGT::SpMatRefType& spmat,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        // OutSource the printer
        ox << prt_prefix << std::endl
           << "->SIZE : [" << spmat.n_rows << " x " << spmat.n_cols << " : "
           << spmat.n_elem << "]" << std::endl;
        brief_print<typename EGT::SpMatRefType::IMatType,
                    typename EGT::IMatType>(spmat.coo, "->COO", ox);
        brief_print<typename EGT::SpMatRefType::IMatType,
                    typename EGT::IMatType>(spmat.csrptr, "->CSRPTR", ox);
        brief_print(spmat.data, "->DATA", ox);
    }

    static void print(const typename EGT::SparseCfgType& spcfg,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        //  OutSource the printer
        ox << prt_prefix << std::endl
           << "->NCONFIG : " << spcfg.n_config << std::endl;
        brief_print<typename EGT::SparseCfgType::CVecType,
                    typename EGT::IVecType>(spcfg.counts, "->CTS:", ox);
        brief_print(spcfg.unique, "->UNQ:", ox);
        brief_print<typename EGT::SparseCfgType::LVecType,
                    typename EGT::IVecType>(spcfg.lookup, "->LKP:", ox);
    }

    static void print(const typename EIMT::model_params_t& mpx,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        ox << prt_prefix << std::endl;
        brief_print(mpx._gamma, "->gamma", ox);
        brief_print(mpx._delta, "->delta", ox);
    }

    static void print(const typename EIMT::batch_adjust_t& bax,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        ox << prt_prefix << std::endl;
        brief_print(bax.denom, "->denom", ox);
        brief_print(bax.numer, "->numer", ox);
    }

    static void print(const EIMT& egimpl, std::ostream& ox = std::cout) {
        ox << "NTHREADS : " << egimpl.nthreads() << std::endl
           << "NGENES   : " << egimpl.nfeatures() << std::endl
           << "NCELLS   : " << egimpl.ncells() << std::endl
           << "NCOVAR   : " << egimpl.ncovariates() << std::endl
           << "NBATCH   : " << egimpl.nbatches() << std::endl
           << "NCONFG   : " << egimpl.nconfig() << std::endl
           << "NEBCFG   : " << egimpl.nebtcfg() << std::endl
           << "NARRAY   : " << egimpl.narray() << std::endl;
        print(egimpl.adata(), "ADATA:", ox);
        print(egimpl.config(), "CFG:", ox);
        print(egimpl.ebt_config(), "EBTCFG:", ox);
        brief_print(egimpl.design(), "DSGN :", ox);
        brief_print(egimpl.batch_counts(), "BCTX:", ox);
        brief_print(egimpl.batch_design(), "EBDSGN :", ox);
        brief_print(egimpl.ebt_counts(), "EBCTX:", ox);
    }
};

template<typename CTMatType=ScementEigenInterface::FMatType,
         bool transpose_output=false>
class ScementEigen : public ScementEigenInterface {
  public:
    //
    using CellSelectType = Eigen::Matrix<uint8_t, 1, Eigen::Dynamic>;
    using FeatSelectType = Eigen::Matrix<uint8_t, Eigen::Dynamic, 1>;
    using SpMatType = COOCSR<FMatType, IMatType, ElemType, IntType>;
    using SpMatRefType = COOCSRRef<SpMatType>;
    using SparseCfgType = SparseCfg<IMatType, IVecType, IVecType>;
    using SparseCfgSliceType = SparseCfgSlice<SparseCfgType, CellSelectType>;
    //
    using ScementEigenImpl = ScementImpl<ScementEigen, CTMatType, transpose_output>;
    using ScementEigenPrtr =
        ScementEigenPrintHelper<ScementEigen, ScementEigenImpl>;
    using BECTMatType = CTMatType;

  protected:
    SpMatType _adata;
    SparseCfgType _config;
    SparseCfgType _ebt_config;
    FMatType _batch_counts;
    FMatType _design;
    FMatType _batch_design;
    ScementEigenImpl _impl;
    ScementEigenPrtr _prn;

  public:
    ScementEigenImpl& impl() { return _impl; }
    ScementEigenPrtr& prn() { return _prn; }

    explicit ScementEigen(DataInterface<ScementEigen>* ploader)
        : _adata(ploader->adata()), _config(ploader->config()),
          _ebt_config(ploader->ebt_config()),
          _batch_counts(ploader->batch_counts()), _design(ploader->design()),
          _batch_design(_design.leftCols(_batch_counts.size())),
          _impl(_adata, _design, _batch_design, _config, _ebt_config,
                _batch_counts, _design.cols(), _batch_counts.size()) {}

    void print(std::ostream& ox = std::cout) const { _prn.print(_impl, ox); }

    template <typename PRDT>
    void print(const PRDT& prx, const std::string header,
               std::ostream& ox = std::cout) const {
        _prn.print(prx, header, ox);
    }
};

template<typename CTMatType=ScementEigenInterface::FMatType,
         bool transpose_output=false>
class ScementEigenMap : public ScementEigenInterface {
  public:
    //
    using CellSelectType = Eigen::Matrix<uint8_t, 1, Eigen::Dynamic>;
    using FeatSelectType = Eigen::Matrix<uint8_t, Eigen::Dynamic, 1>;
    //
    using SpMatType = COOCSR<FMatType, MapIMatType, ElemType, IntType>;
    using SparseCfgType = SparseCfg<IMatType, IVecType, MapIVecType>;
    using SpMatRefType = COOCSRRef<SpMatType>;
    using SparseCfgSliceType = SparseCfgSlice<SparseCfgType, CellSelectType>;
    //
    using ScementEigMapImpl = ScementImpl<ScementEigenMap, CTMatType, transpose_output>;
    using ScementEigenMapPrtr =
        ScementEigenPrintHelper<ScementEigenMap, ScementEigMapImpl>;
    using BECTMatType = CTMatType;

  protected:
    IndexType _nfeatures, _ncells, _nelems;
    SpMatType _adata;
    SparseCfgType _config;
    SparseCfgType _ebt_config;
    FMatType _batch_counts;
    FMatType _design;
    FMatType _batch_design;

    ScementEigMapImpl _impl;
    ScementEigenMapPrtr _prn;
  public:
    ScementEigMapImpl& impl() { return _impl; }
    ScementEigenMapPrtr& prn() { return _prn; } 

    explicit ScementEigenMap(DataInterface<ScementEigenMap>* ploader)
        : _adata(ploader->adata()), _config(ploader->config()),
          _ebt_config(ploader->ebt_config()),
          _batch_counts(ploader->batch_counts()), _design(ploader->design()),
          _batch_design(_design.leftCols(_batch_counts.size())),
          _impl(_adata, _design, _batch_design, _config, _ebt_config,
                _batch_counts, _design.cols(), _batch_counts.size()) {}

    void print(std::ostream& ox = std::cout) const {
        _prn.print(_impl, ox);
    }
    template <typename PRDT>
    void print(const PRDT& prx, const std::string header,
               std::ostream& ox = std::cout) const {
        _prn.print(prx, header, ox);
    }
};

#endif  // !SCEMENT_ARMA_HPP
