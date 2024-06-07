//
// Copyright [2024]
//
#ifndef SCEMENT_IMPL_ARMA_HPP
#define SCEMENT_IMPL_ARMA_HPP

#include <cstdint>
#include <string>

#include "scement/coocsr.hpp"
#include "scement/data_if.hpp"
#include "scement/impl_base.hpp"
#include <armadillo>

class ScementArmaInterface {
public:
    using ElemType = float;
    using IntType = int32_t;
    using UIntType = uint32_t;
    using IndexType = uint32_t;
    using SelectType = uint8_t;
    // 
    using FMatType = arma::fmat;
    using FVecType = arma::fvec;
    using IVecType = arma::Col<int32_t>;
    using FRowType = arma::Row<float>;
    using IMatType = arma::Mat<int32_t>;
    using UMatType = arma::Mat<uint32_t>;
    // 
    static inline ElemType pinv_tolerance() { return 1e-8; }
    static inline FMatType transpose(const FMatType& infm) { return infm.t(); }
    static inline IMatType transpose(const IMatType& infm) { return infm.t(); }
    static inline FMatType imat2fmat(const IMatType& infm) {
        return arma::conv_to<FMatType>::from(infm);
    }
    static inline FRowType colsum(const FMatType& infm) {
        return arma::sum(infm, 0);
    }
    static inline FVecType colmean(const FMatType& infm) {
        return arma::mean(infm, 0).as_col();
    }
    static inline FMatType diag(const FRowType& rx) {
        return arma::diagmat(rx);
    }
    static inline FMatType pseudo_inverse(const FMatType& infm) {
        return arma::pinv(infm, pinv_tolerance());
    }
    static inline FMatType head_columns(const FMatType& infm,
                                        std::size_t ncols) {
        return infm.head_cols(ncols);
    }
    static inline FVecType colvar(const FMatType& infm) {
        return arma::var(infm, 0).as_col();
    }
    static inline ElemType mean(const FMatType& infm) {
        return arma::mean(infm.as_col());
    }
    static inline ElemType var(const FMatType& infm) {
        return arma::var(infm.as_col());
    }
    static inline IntType sum(const IMatType& infm) { return arma::accu(infm); }
    //
    //
    static inline FMatType gram_matrix(const FMatType& infm) {
        return infm.t() * infm;
    }
    static inline FMatType bhat(const FMatType& des_adt,
                                const FMatType& des_gram) {
        // FMatType des_pv = pseudo_inverse(des_gram);  // NCVR x NCVR
        // _B_hat = des_adt * transpose(des_pv);        // NFEATURESxNCVR
        return des_adt * arma::pinv(des_gram, pinv_tolerance()).t();
    }
    static inline FMatType cfg_bhat(const FMatType& bhx,
                                    const IMatType& cfg_unq) {
        // _B_hat *
        //  ScementData::imat2fmat(ScementData::transpose(_config.unique));
        return bhx * cfg_unq.t();
    }
    static inline FMatType grand_mean(const FMatType& bhx,
                                      const FMatType& ctx) {
        // _grand_mean =
        //     ScementData::head_columns(_B_hat, _nbatches) *
        //   (ScementData::transpose(_batch_counts) / _narray);  // NFEATURESx1
        return bhx.head_cols(ctx.size()) * (ctx / arma::accu(ctx));
    }

    static inline FMatType ghat_base(const FMatType& bdes) {
        return bdes *
               arma::pinv(arma::diagmat(arma::sum(bdes, 0)), pinv_tolerance());
    }
    static inline FMatType ghat_mprime(const FMatType& mp,
                                       const IMatType& ctx) {
        return mp * ctx;
        //_m_prime * SDT::imat2fmat(_ebt_counts);  // NFEATURESxNBATCH
    }
};

template<typename AMT, typename AMIMT>
class ScementArmaPrintHelper {
public:
    static void print(const typename AMT::IMatType& fmat,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        fmat.brief_print(ox, prt_prefix);
    }

    static void print(const typename AMT::FMatType& fmat,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        fmat.brief_print(ox, prt_prefix);
    }

    static void print(const typename AMT::FVecType& xmat,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        xmat.brief_print(ox, prt_prefix);
    }

   static void print(const typename AMT::SpMatType& spmat,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        ox << prt_prefix << std::endl 
           << "->SIZE : [" << spmat.n_rows << " x " << spmat.n_cols 
           << " : " << spmat.n_elem << "]" << std::endl;
        spmat.coo.brief_print(ox, "->COO");
        spmat.csrptr.brief_print(ox, "->CSRPTR");
        spmat.data.brief_print(ox, "->DATA");
    }

    static void print(const typename AMT::SpMatRefType& spmat,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout)  {
        ox << prt_prefix << std::endl 
           << "->SIZE : [" << spmat.n_rows << " x " << spmat.n_cols 
           << " : " << spmat.n_elem << "]" << std::endl;
        spmat.coo.brief_print(ox, "->COO");
        spmat.csrptr.brief_print(ox, "->CSRPTR");
        spmat.data.brief_print(ox, "->DATA");
    }

    static void print(const typename AMT::SparseCfgType& spcfg,
                      const std::string prt_prefix, 
                      std::ostream& ox = std::cout)  {
        ox << prt_prefix << std::endl
           << "->NCONFIG : " << spcfg.n_config << std::endl;
        spcfg.counts.brief_print(ox, "->CTS:");
        spcfg.unique.brief_print(ox, "->UNQ:");
        spcfg.lookup.brief_print(ox, "->LKP:");
    }

    static void print(const typename AMIMT::model_params_t& mpx,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout)  {
        ox << prt_prefix << std::endl;
        mpx._gamma.brief_print(ox, "->gamma");
        mpx._delta.brief_print(ox, "->delta");
    }

    static void print(const typename AMIMT::batch_adjust_t& bax,
                      const std::string prt_prefix,
                      std::ostream& ox = std::cout) {
        ox << prt_prefix << std::endl;
        bax.denom.brief_print(ox, "->denom");
        bax.numer.brief_print(ox, "->numer");
    }

    static void print(const AMIMT& arimpl, std::ostream& ox = std::cout) {
        ox << "NTHREADS : " << arimpl.nthreads() << std::endl
           << "NGENES   : " << arimpl.nfeatures() << std::endl
           << "NCELLS   : " << arimpl.ncells() << std::endl
           << "NCOVAR   : " << arimpl.ncovariates() << std::endl
           << "NBATCH   : " << arimpl.nbatches() << std::endl
           << "NCONFG   : " << arimpl.nconfig() << std::endl
           << "NEBCFG   : " << arimpl.nebtcfg() << std::endl
           << "NARRAY   : " << arimpl.narray() << std::endl;
        print(arimpl.adata(), "ADATA:", ox);
        print(arimpl.config(), "CFG:", ox);
        print(arimpl.ebt_config(), "EBTCFG:", ox);
        arimpl.design().brief_print(ox, "DSGN :");
        arimpl.batch_counts().brief_print(ox, "BCTX:");
        arimpl.batch_design().brief_print(ox, "EBDSGN :");
        arimpl.ebt_counts().brief_print(ox, "EBCTX:");
    }
};

template<typename CTMatType=ScementArmaInterface::FMatType,
         bool transpose_output=false>
class ScementArma : public ScementArmaInterface {
  public:
    // part of interface
    using CellSelectType = arma::Row<uint8_t>;
    using FeatSelectType = arma::Col<uint8_t>;
    //
    using SpMatType = COOCSR<FMatType, IMatType, ElemType, IntType>;
    using SpMatRefType = COOCSRRef<SpMatType>;
    using SparseCfgType = SparseCfg<IMatType, IVecType, IVecType>;
    using SparseCfgSliceType = SparseCfgSlice<SparseCfgType, CellSelectType>;
    //
    using ScementArmaImpl = ScementImpl<ScementArma, CTMatType, transpose_output>;
    using ScementArmaPrtr = ScementArmaPrintHelper<ScementArma, ScementArmaImpl>;
    using BECTMatType = CTMatType;
  private:
    SpMatType _adata;
    SparseCfgType _config;
    SparseCfgType _ebt_config;
    FMatType _batch_counts;
    FMatType _design;
    FMatType _batch_design;
    ScementArmaImpl _impl;
    ScementArmaPrtr _prn;
  public:
    ScementArmaImpl& impl() { return _impl; }
    ScementArmaPrtr& prn() { return _prn; }

    explicit ScementArma(DataInterface<ScementArma>* ploader)
        : _adata(ploader->adata()), _config(ploader->config()),
          _ebt_config(ploader->ebt_config()),
          _batch_counts(ploader->batch_counts()), _design(ploader->design()),
          _batch_design(_design.head_cols(_batch_counts.size())),
          _impl(_adata, _design, _batch_design, _config, _ebt_config,
                _batch_counts, _design.n_cols, _batch_counts.size()) {}

    void print(std::ostream& ox = std::cout) const {
        _prn.print(_impl, ox);
    }
    template <typename PRDT>
    void print(const PRDT& prx, const std::string header,
               std::ostream& ox = std::cout) const {
        _prn.print(prx, header, ox);
    }
};

#endif  // !SCEMENT_IMPL_ARMA_HPP
