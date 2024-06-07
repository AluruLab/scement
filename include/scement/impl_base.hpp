#ifndef IMPL_BASE_HPP
#define IMPL_BASE_HPP

#include "scement/coocsr.hpp"
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

template <typename SDT, typename CTMatType = typename SDT::FMatType,
          bool transposed_output = false>
class ScementImpl {
  public:
    using ElemType = typename SDT::ElemType;
    using IndexType = typename SDT::IndexType;
    using IntType = typename SDT::IntType;
    //
    using FMatType = typename SDT::FMatType;
    using IMatType = typename SDT::IMatType;
    using FVecType = typename SDT::FVecType;
    using IVecType = typename SDT::IVecType;
    //
    using SelectType = typename SDT::SelectType;
    using CellSelectType = typename SDT::CellSelectType;
    using FeatSelectType = typename SDT::FeatSelectType;
    using SpMatType = typename SDT::SpMatType;
    using SpMatRefType = typename SDT::SpMatRefType;
    using SparseCfgType = typename SDT::SparseCfgType;
    using SparseCfgSliceType = typename SDT::SparseCfgSliceType;
    //
    struct model_params_t {
        FMatType _gamma;  // NFEATURESx1
        FMatType _delta;  // NFEATURESx1
    };
    struct batch_adjust_t {
        FMatType numer;  // NFEATURESxNEBCFG (same as _mprime)
        FVecType denom;  // NFEATURESx1
    };
    //
    const ElemType pinv_tolerance = 1e-8;

  private:
    int _nthreads;
    // Inputs
    const SpMatType& _adata;        // NFEATURESxNCELLS
    const FMatType& _design;        // NCELLSxNCOVRS
    const FMatType& _batch_design;  // NCELLSxNBATCH
    const SparseCfgType& _config;   // C:NCVCFGx1; U:NCVCFGxNCOVARS; L:NCELLSx1
    const SparseCfgType&
        _ebt_config;         // C:NEBCFGx1; U:NEBCFGxNCOVARS; L:NCELLSx1
    FMatType _batch_counts;  // NBATCHx1
    IndexType _nfeatures, _ncells;
    IndexType _ncovariates, _nbatches, _nconfig, _nebtcfg;
    // Computed Data
    ElemType _narray;                  // sum()
    FVecType _var_pooled;              // NFEATURESx1
    FMatType _cfgB_hat;                // NFEATURESxNCOVARS
    FMatType _B_hat;                   // NFEATURESxNCOVARS
    FMatType _standB_hat;              // NEBCFGxNFEATURES
    FMatType _grand_mean;              // NFEATURESxNEBCFG
    FMatType _stand_mean;              // NFEATURESxNEBCFG
    SpMatRefType _d_prime;             // NFEATURESxNCELLS
    FMatType _m_prime;                 // NFEATURESxNEBCFG
    FMatType _ghat_pfx;                // NCELLSxNBATCH
    FMatType _ghat_dprime;             // NFEATURESxNBATCH
    FMatType _ghat_mprime;             // NFETURESxNBATCH
    FMatType _ghat;                    // NFETURESxNBATCH
    IMatType _ebt_counts;              // NEBCFGxNBATCH
    std::vector<FMatType> _delta_hat;  // NBATCH[NFEATURESx1]
    FVecType _gamma_bar;               // NBATCHx1
    FVecType _t2;                      // NBATCHx1
    FVecType _aprior;                  // NBATCHx1
    FVecType _bprior;                  // NBATCHx1
    //
    CellSelectType _cell_select;      // NCELLSx1
    CellSelectType _cell_select_all;  // NCELLSx1
    FeatSelectType _feat_select;      // NFEATURESx1
    FeatSelectType _feat_select_all;  // NFEATURESx1
    //
    std::vector<model_params_t> _param_star;
    std::vector<batch_adjust_t> _batch_adj;

  public:
    // Getter functions
    inline int nthreads() const { return _nthreads; }
    inline IndexType nfeatures() const { return _nfeatures; }
    inline IndexType ncells() const { return _ncells; }
    inline IndexType ncovariates() const { return _ncovariates; }
    inline IndexType nbatches() const { return _nbatches; }
    inline IndexType nebtcfg() const { return _nebtcfg; }
    inline IndexType nconfig() const { return _nconfig; }
    inline ElemType narray() const { return _narray; }
    inline const SpMatType& adata() const { return _adata; }
    inline const SparseCfgType& config() const { return _config; }
    inline const SparseCfgType& ebt_config() const { return _ebt_config; }
    inline const FMatType& batch_counts() const { return _batch_counts; }
    inline const FMatType& design() const { return _design; }
    inline const FMatType& batch_design() const { return _batch_design; }
    inline const FVecType& var_pooled() const { return _var_pooled; }
    inline const FMatType& cfgB_hat() const { return _cfgB_hat; }
    inline const FMatType& B_hat() const { return _B_hat; }
    inline const FMatType& standB_hat() const { return _standB_hat; }
    inline const FMatType& grand_mean() const { return _grand_mean; }
    inline const FMatType& std_mean() const { return _stand_mean; }
    inline const SpMatRefType& d_prime() const { return _d_prime; }
    inline const FMatType& m_prime() const { return _m_prime; }
    inline const FMatType& ghat_pfx() const { return _ghat_pfx; }
    inline const FMatType& ghat_dprime() const { return _ghat_dprime; }
    inline const FMatType& ghat_mprime() const { return _ghat_mprime; }
    inline const FMatType& ghat() const { return _ghat; }
    inline const IMatType& ebt_counts() const { return _ebt_counts; }
    inline const std::vector<FMatType>& delta_hat() const { return _delta_hat; }
    inline const FVecType& gamma_bar() const { return _gamma_bar; }
    inline const FVecType& t2() const { return _t2; }
    inline const FVecType& aprior() const { return _aprior; }
    inline const FVecType& bprior() const { return _bprior; }
    inline const std::vector<model_params_t>& param_star() const {
        return _param_star;
    }
    inline const std::vector<batch_adjust_t>& batch_adj() const {
        return _batch_adj;
    }
    //

  private:
    // Data Initialization
    void init() {
        // _nfeatures = _adata.n_rows;
        // _ncells = _adata.n_cols;
        // _nebtcfg = _ebt_config.n_config;
        // Selection Flag Array
        _cell_select = CellSelectType(_ncells);
        _feat_select = FeatSelectType(_nfeatures);
        _cell_select_all = CellSelectType(_ncells);
        _feat_select_all = FeatSelectType(_nfeatures);
        //
#pragma omp parallel for
        for (IndexType i = 0; i < _ncells; i++) {
            _cell_select(i) = 0;
            _cell_select_all(i) = 1;
        }
#pragma omp parallel for
        for (IndexType i = 0; i < _nfeatures; i++) {
            _feat_select(i) = 0;
            _feat_select_all(i) = 1;
        }
        // Num of total entries
        _narray = 0.0;
        for (IndexType i = 0; i < _batch_counts.size(); i++) {
            _narray += _batch_counts(i);
        }
        // Param Initializations
        for (IndexType bx = 0; bx < _nbatches; bx++) {
            _aprior(bx) = _bprior(bx) = 0;
        }
        // Config
        for (IndexType i = 0; i < _ebt_counts.size(); i++) {
            _ebt_counts(i) = 0;
        }
        // Batch Design
        for (IndexType jx = 0; jx < _nebtcfg; jx++) {
            for (IndexType bx = 0; bx < _nbatches; bx++) {
                ElemType sum = 0.0;
#pragma omp parallel for reduction(+ : sum)
                for (auto cx = 0; cx < _ncells; cx++) {
                    ElemType ind_val =
                        static_cast<ElemType>(_ebt_config.lookup(cx) == jx);
                    sum += ind_val * _batch_design(cx, bx);
                }
                _ebt_counts(jx, bx) = sum;
            }
        }
    }

  private:
    // 0. B_hat and cfgB_hat
    void bhat() {
        //
        FMatType des_gram = SDT::gram_matrix(_design);  // NCVR x NCVR
        //
        FMatType des_adt = mult_coocsr_dense<FMatType, ElemType>(
            _adata, _design, _adata.n_rows, _ncovariates);
        // FMatType des_pv =
        //     SDT::pseudo_inverse(udes_gram);  // NCVR x NCVR
        // _B_hat = des_adt * SDT::transpose(des_pv);  // NFEATURESxNCVR
        _B_hat = SDT::bhat(des_adt, des_gram);              // NFEATURESxNCVR
        _cfgB_hat = SDT::cfg_bhat(_B_hat, _config.unique);  // NFEATURESxNCVR
        //     _B_hat *
        //    SDT::imat2fmat(SDT::transpose(_config.unique));
    }

    // 1 pooled variance
    void compute_var_pooled() {
        _var_pooled = row_sum_of_squares<FVecType>(
            _adata, _cfgB_hat, SparseCfgSliceType(_config, _cell_select_all));
        _var_pooled /= _narray;  // NFEATURESx1
    }

    // 2 standardized mean
    void stand_mean() {
        _grand_mean =
            SDT::grand_mean(_B_hat, _batch_counts);  // NFEATURESxNEBCFG
        _standB_hat =
            SDT::cfg_bhat(_B_hat, _ebt_config.unique);  // NFEATURESxNEBCFG
        //
        _stand_mean = FMatType(_nfeatures, _nebtcfg);
#pragma omp parallel for
        for (IndexType col = 0; col < _nebtcfg; col++) {
#pragma omp parallel for
            for (IndexType row = 0; row < _nfeatures; row++) {
                _stand_mean(row, col) =
                    _standB_hat(row, col) + _grand_mean(row);
            }
        }
        // _stand_mean = _standB_hat.each_col() + _grand_mean;  //
        // NFEATURESxNEBCFG
    }

    // 3 dprime and mprime
    void standardize_data() {
        FVecType col_avp(_nfeatures);  // NFEATURESx1
        // col_avp = 1.0 / SDT::sqrt(_var_pooled);
#pragma omp parallel for
        for (IndexType row = 0; row < _nfeatures; row++) {
            col_avp(row) = 1.0 / std::sqrt(_var_pooled(row));
        }
        //
        elt_mult_coocsr_dense<SpMatType, FVecType, SpMatRefType>(
            _adata, col_avp, &_d_prime);
        //_m_prime = _stand_mean.each_col() % col_avp;  // NFEATURESxNEBCFG
        _m_prime = FMatType(_nfeatures, _nebtcfg);  // NFEATURESxNEBCFG
#pragma omp parallel for
        for (IndexType col = 0; col < _nebtcfg; col++) {
#pragma omp parallel for
            for (IndexType row = 0; row < _nfeatures; row++) {
                _m_prime(row, col) = _stand_mean(row, col) * col_avp(row);
            }
        }

        // # Need to be a bit careful with the zero variance features
        // # just set the zero variance features to zero in the standardized
        // data if np.any(var_pooled == 0):
        //     s_dprime[var_pooled == 0, :] = 0
        //     s_mprime[var_pooled == 0, :] = 0
#pragma omp parallel for
        for (IndexType i = 0; i < _nfeatures; i++) {
            _feat_select(i) = static_cast<int8_t>(_var_pooled(i) == 0);
        }
        _d_prime.set(_feat_select, _cell_select_all, 0.0);
        //
#pragma omp parallel for
        for (IndexType col = 0; col < _nebtcfg; col++) {
#pragma omp parallel for
            for (IndexType row = 0; row < _nfeatures; row++) {
                if (_feat_select(row)) {
                    _m_prime(row, col) = 0.0;
                }
            }
        }
    }

    // 4 gamma hat
    void gamma_hat_prime() {
        // gamma_hat_pfx = np.linalg.pinv(np.diag(
        //     np.array(batch_design.T.sum(axis=1),
        //              dtype=run_dtype))) @ batch_design.T
        _ghat_pfx = SDT::ghat_base(_batch_design);  // NCELLSxNBATCH

        //    _batch_design * SDT::pseudo_inverse(SDT::diag(
        //                        SDT::colsum(_batch_design)));

        // gamma_hat_dprime = (s_dprime @ sp_gamma_hat_pfx.T).A.T
        _ghat_dprime = mult_coocsr_dense<FMatType, float>(
            _d_prime, _ghat_pfx, _d_prime.n_rows,
            _nbatches);  // NFEATURESxNBATCH

        // gamma_hat_mprime = (s_mprime @ batchxstd_counts) / n_batches

        _ghat_mprime =
            SDT::ghat_mprime(_m_prime, _ebt_counts);  // NFEATURESxNBATCH
        // _m_prime * SDT::imat2fmat(_ebt_counts);  // NFEATURESxNBATCH
        //_ghat_mprime.each_row() /= _batch_counts;
#pragma omp parallel for
        for (IndexType col = 0; col < _nbatches; col++) {
#pragma omp parallel for
            for (IndexType row = 0; row < _nfeatures; row++) {
                _ghat_mprime(row, col) /= _batch_counts(col);
            }
        }
        // gamma_hat = gamma_hat_dprime - gamma_hat_mprime.T
        _ghat = _ghat_dprime - _ghat_mprime;  // NFEATURESxNBATCH
    }

    // 5 delta_hat loop
    void multplicative_effect() {
        // delta_hat = []
        // for i, batch_idxs in enumerate(batch_idxvals):
        //     # delta_b = np.zeros() TODO:
        //     cfg_bx = stdcfg_lookup[np.array(batch_design.iloc[:, i],
        //     dtype=bool)]
        //
        //     batch_s_dprime = s_dprime[:, batch_idxs]
        //
        //     mean_bx =
        //     s_mprime + row_mean_sparse(batch_s_dprime,
        //                                          s_mprime.T,
        //                                          batchxstd_counts[:, i],
        //                                          cfg_bx)
        //     var_bx = row_sum_of_squares_sparse(batch_s_dprime,
        //                                        mean_bx,
        //                                        batchxstd_counts[:, i],
        //                                        cfg_bx)
        //     #
        //     delta_hat.append(var_bx / (n_batches[i] - 1))
        // _ebt_counts.print("EBTC:");

        FMatType mean_bx(_nfeatures, _nebtcfg);  // NFEATURESxNEBCFG
        for (auto bx = 0; bx < _nbatches; bx++) {
#pragma omp parallel for
            for (auto cx = 0; cx < _ncells; cx++) {
                _cell_select(cx) =
                    static_cast<SelectType>(_batch_design(cx, bx));
            }
            typename SparseCfgSliceType::CVecType bctx(
                _ebt_counts.col(bx));  // NEBCFGx1
            // bctx.print("EBT CX");
            SparseCfgSliceType batch_cfg(_ebt_config.unique, bctx,
                                         _ebt_config.lookup,
                                         _ebt_config.n_config, _cell_select);
            FVecType rsum =
                coocsr_row_mean<FVecType>(_d_prime, _m_prime, batch_cfg);

            // arma::fmat mean_bx(_m_prime);  // NFEATURESxNEBCFG
#pragma omp parallel for
            for (auto row = 0; row < _nfeatures; row++) {
                for (auto col = 0; col < _nebtcfg; col++) {
                    mean_bx(row, col) = _m_prime(row, col) + rsum(row);
                }
            }

            // mean_bx.each_col() += rsum;
            _delta_hat[bx] =  // NFEATURESx1
                row_sum_of_squares<FVecType>(_d_prime, mean_bx, batch_cfg);
#pragma omp parallel for
            for (auto fx = 0; fx < _nfeatures; fx++) {
                _delta_hat[bx](fx) /= _batch_counts(bx);
            }
            // var_bx /= _batch_counts[bx];  // NFEATURESx1
        }
    }

    // 5 fix hyperparms
    void fix_hyperparams() {
        //
        // gamma_bar = gamma_hat.mean(axis=1)
        _gamma_bar = SDT::colmean(_ghat);
        // t2 = gamma_hat.var(axis=1)
        _t2 = SDT::colvar(_ghat);
        //
        // # a_prior and b_prior are the priors on lambda and theta from Johnson
        // and Li (2006)
        // NCOVARS
        // m = delta_hat.mean()
        // s2 = delta_hat.var()
        // _aprior: (2 * s2 + m**2) / s2
        // _bprior: (m * s2 + m**3) / s2
        for (auto bx = 0; bx < _nbatches; bx++) {
            float m = SDT::mean(_delta_hat[bx]);
            float s2 = SDT::var(_delta_hat[bx]);
            _aprior(bx) = (2 * s2 + (m * m)) / s2;
            _bprior(bx) = (m * s2 + (m * m * m)) / s2;
            if (std::isinf(_bprior(bx)))
                _bprior(bx) = 2 * _aprior(bx);
        }
    }

    model_params_t _it_sol(const SparseCfgSliceType& sub_config,
                           const model_params_t& e_hat, int32_t bt_counts,
                           float g_bar, float t2, float a, float b,
                           float conv = 0.0001) {

        // g_old = g_hat.copy()
        // d_old = d_hat.copy()
        // n = np.ones((s_dprime.shape[0],)) * s_dprime.shape[1]
        unsigned count = 0;
        FVecType n(_d_prime.n_rows);
        for (IndexType i = 0; i < _d_prime.n_rows; i++) {
            n(i) = 1;
        }
        n *= static_cast<float>(bt_counts);
        model_params_t e_old = e_hat;
        model_params_t e_new = e_old;
        // arma::fvec change(1, arma::fill::ones);
        float change = 1.0;
        //
        //_m_prime = FMatType(_nfeatures, _nebtcfg);  // NFEATURESxNEBCFG
        FMatType g_new_mp(_nfeatures, _nebtcfg);
        while (change > conv) {
            // g_new = (t2 * n * g_hat + d_old * g_bar) / (t2 * n + d_old)
#pragma omp parallel for
            for (IndexType i = 0; i < e_new._gamma.size(); i++) {
                e_new._gamma(i) =
                    (e_hat._gamma(i) * n(i) * t2) + (e_old._delta(i) * g_bar);
                e_new._gamma(i) /= e_old._delta(i) + (n(i) * t2);
            }
            // --- CPP ---
            // e_new._gamma = ((e_hat._gamma % n) * t2) + (e_old._delta *
            // g_bar); e_new._gamma /= e_old._delta + (n * t2);
            // --- CPP ---
            // g_new_mp = s_mprime - g_new.reshape((g_new.shape[0], 1))
            // --- CPP ---
            // arma::fmat g_new_mp = _m_prime.each_col() - e_new._gamma;
            // --- CPP ---
            // Same dimensionsa _m_prime (_nfeatures, _nebtcfg)
            for (IndexType col = 0; col < _nebtcfg; col++) {
#pragma omp parallel for
                for (IndexType row = 0; row < _nfeatures; row++) {
                    g_new_mp(row, col) = _m_prime(row, col) - e_new._gamma(row);
                }
            }
            //
            // sum2 = row_sum_of_squares_sparse(s_dprime, g_new_mp,
            //                                  stdcfg_counts, stdcfg_lookup)
            FVecType sum2 =
                row_sum_of_squares<FVecType>(_d_prime, g_new_mp, sub_config);
            // sum2 = sum2.flatten()
            // arma::fvec sum2 = rsum2.as_col();
            // d_new = (0.5 * sum2 + b) / (n / 2.0 + a - 1.0)
            // --- CPP ---
            // e_new._delta = (0.5 * sum2 + b);
            // e_new._delta /= (n / 2.0) + a - 1.0;
            // --- CPP ---
#pragma omp parallel for
            for (auto i = 0; i < e_new._delta.size(); i++) {
                e_new._delta(i) = 0.5 * sum2(i) + b;
                e_new._delta(i) /= (n(i) / 2.0) + a - 1.0;
            }
            // gsel = (g_old != 0)
            // dsel = (d_old != 0)
            // change = max((abs(g_new[gsel] - g_old[gsel]) /
            // g_old[gsel]).max(),
            //         (abs(d_new[dsel] - d_old[dsel]) / d_old[dsel]).max())
            float gumax = 0.0, dumax = 0.0;
#pragma omp parallel for reduction(max : gumax)
            for (auto i = 0; i < e_old._gamma.size(); i++) {
                float gupd = 0.0;
                if (e_old._gamma(i) != 0) {
                    gupd = std::abs(e_new._gamma(i) - e_old._gamma(i)) /
                           e_old._gamma(i);
                }
                gumax = gumax > gupd ? gumax : gupd;
            }
#pragma omp parallel for reduction(max : dumax)
            for (auto i = 0; i < e_old._delta.size(); i++) {
                float dupd = 0.0;
                if (e_old._delta(i) != 0) {
                    dupd = std::abs(e_new._delta(i) - e_old._delta(i)) /
                           e_old._delta(i);
                }
                dumax = dumax > dupd ? dumax : dupd;
            }
            change = std::max(gumax, dumax);
            // change.print("CHG");
            // g_old = g_new  # .copy()
            // d_old = d_new  # .copy()
            e_old = e_new;
            // count = count + 1
            count = count + 1;
        }
        return e_new;
    }

    // 6 Compute parametric adjustments
    void parametric_adjustments() {
        for (IndexType i = 0; i < _nbatches; i++) {
#pragma omp parallel for
            for (IndexType kx = 0; kx < _ncells; kx++) {
                _cell_select(kx) =
                    static_cast<SelectType>(_batch_design(kx, i));
            }
            IVecType bctx(_ebt_counts.col(i));
            IntType bt_counts = SDT::sum(bctx);
            SparseCfgSliceType batch_cfg(_ebt_config.unique, bctx,
                                         _ebt_config.lookup,
                                         _ebt_config.n_config, _cell_select);
            // gamma, delta = _it_sol_sparse(
            //     s_dprime[:, batch_idxs], s_mprime, batchxstd_counts[:, i],
            //     stdcfg_lookup[batch_idxs], gamma_hat[i],
            //     delta_hat[i].flatten(), gamma_bar[i], t2[i],
            //     float(a_prior[i]), float(b_prior[i]),)
            model_params_t ehat;
            ehat._gamma = _ghat.col(i);
            ehat._delta = _delta_hat[i];
            _param_star[i] = _it_sol(batch_cfg, ehat, bt_counts, _gamma_bar(i),
                                     _t2(i), _aprior(i), _bprior(i));
            // gamma_star.append(gamma)
            // delta_star.append(delta)
        }
    }

    void adjust_data() {
        // logg.info("--> Adjusting data")
        // batch_adjustment = []

        // for j, batch_idxs in enumerate(batch_idxvals):
        for (auto i = 0; i < _nbatches; i++) {
            //  dsq = np.sqrt(delta_star[j, :])
            // _batch_adj[i].denom = SDT::sqrt(_param_star[i]._delta);
            _batch_adj[i].denom = _param_star[i]._delta;
#pragma omp parallel for
            for (IndexType row = 0; row < _batch_adj[i].denom.size(); row++) {
                _batch_adj[i].denom(row) =
                    std::sqrt(_param_star[i]._delta(row));
            }
            //  denom = dsq.reshape((len(dsq), 1))
            //  cfg_gamma_star = gamma_star[j, :].reshape((n_features, 1))
            _batch_adj[i].numer = _m_prime;
            // Same dimensionsa _m_prime (_nfeatures, _nebtcfg)

#pragma omp parallel for
            for (IndexType col = 0; col < _nebtcfg; col++) {
#pragma omp parallel for
                for (IndexType row = 0; row < _nfeatures; row++) {

                    _batch_adj[i].numer(row, col) += _param_star[i]._gamma(row);
                }
            }
            // _batch_adj[i].numer.each_col() += _param_star[i]._gamma;
            //  batch_adjustment.append((cfg_gamma_star + s_mprime, denom))
        }
        // log_mem_usage('Compuated Batch Adjustments for all Batches')
    }

    template <bool enabled = transposed_output>
    inline void
    set_value(CTMatType* pbc_matrix,
              IndexType row, IndexType col, ElemType val,
              typename std::enable_if<enabled>::type* = 0) {
        (*pbc_matrix)(col, row) = val;
    }
    template <bool enabled = transposed_output>
    inline void
    set_value(CTMatType* pbc_matrix,
              IndexType row, IndexType col, ElemType val,
              typename std::enable_if<!enabled>::type* = 0) {
        (*pbc_matrix)(row, col) = val;
    }

    template <bool enabled = transposed_output>
    inline ElemType
    value_at(CTMatType* pbc_matrix,
              IndexType row, IndexType col, 
              typename std::enable_if<enabled>::type* = 0) {
        return (*pbc_matrix)(col, row);
    }
    template <bool enabled = transposed_output>
    inline ElemType
    value_at(CTMatType* pbc_matrix,
              IndexType row, IndexType col, 
              typename std::enable_if<!enabled>::type* = 0) {
        return (*pbc_matrix)(row, col);
    }

  public:
    CTMatType* construct_bcdata(CTMatType* pbc_matrix) {
        assert(pbc_matrix->size() == (_nfeatures * _ncells));
        //
        // bayesdata = s_dprime.A
        // RFMatType bc_matrix(_nfeatures, _ncells);
#pragma omp parallel for
        for (IndexType col = 0; col < _ncells; col++) {
#pragma omp parallel for
            for (IndexType row = 0; row < _nfeatures; row++) {
                // (*pbc_matrix)(row, col) = 0;
                set_value(pbc_matrix, row, col, 0.0);
            }
        }
#pragma omp parallel for
        for (IndexType sdx = 0; sdx < _d_prime.n_elem; sdx++) {
            IndexType row_id = _d_prime.coo(sdx, 0);
            IndexType col_id = _d_prime.coo(sdx, 1);
            // (*pbc_matrix)(row_id, col_id) = _d_prime.data(sdx);
            set_value(pbc_matrix, row_id, col_id, _d_prime.data(sdx));
        }

        // log_mem_usage('Covert to data to Dense')
        // for i in range(n_stdcfgs):
        // for j, batch_idxs in enumerate(batch_idxvals):
        // cell_sel = (stdcfg_lookup[batch_idxs] == i)
        // numer_adj, denom = batch_adjustment[j]
        // if np.any(cell_sel):
        //     bayesdata[:, batch_idxs[cell_sel]] -= numer_adj[:,
        //     i].reshape((
        //         n_features, 1))
        //     bayesdata[:, batch_idxs[cell_sel]] /= denom

        // vpsq = np.sqrt(var_pooled).reshape((len(var_pooled), 1))
        FVecType vpsq(_nfeatures);  // SDT::sqrt(_var_pooled);
#pragma omp parallel for
        for (IndexType row = 0; row < _nfeatures; row++) {
            vpsq(row) = std::sqrt(_var_pooled(row));
        }

#pragma omp parallel for
        for (auto cell_id = 0; cell_id < _ncells; cell_id++) {
            IntType batch_id = -1;
            for (IntType bid = 0; bid < _nbatches; bid++) {
                if (_batch_design(cell_id, bid) == 1) {
                    batch_id = bid;
                    break;
                }
            }
            if (batch_id == -1) {
                continue;
            }
            IntType ex = _ebt_config.lookup(cell_id);

#pragma omp parallel for
            for (auto feat_id = 0; feat_id < _nfeatures; feat_id++) {
                // (*pbc_matrix)(feat_id, cell_id) -=
                //     _batch_adj[batch_id].numer(feat_id, ex);
                // (*pbc_matrix)(feat_id, cell_id) /=
                //     _batch_adj[batch_id].denom(feat_id);
                set_value(pbc_matrix, feat_id, cell_id,
                    value_at(pbc_matrix, feat_id, cell_id) -
                    _batch_adj[batch_id].numer(feat_id, ex));
                set_value(pbc_matrix, feat_id, cell_id,
                    value_at(pbc_matrix, feat_id, cell_id) /
                    _batch_adj[batch_id].denom(feat_id));
            }
            // _bcdata.col(cell_id) -= _batch_adj[batch_id].numer.col(ex);
            // _bcdata.col(cell_id) /= _batch_adj[batch_id].denom;
        }
        // _bcdata.brief_print();
        // bayesdata = np.multiply(bayesdata, vpsq)
#pragma omp parallel for
        for (auto cell_id = 0; cell_id < _ncells; cell_id++) {
#pragma omp parallel for
            for (auto feat_id = 0; feat_id < _nfeatures; feat_id++) {
                // (*pbc_matrix)(feat_id, cell_id) *= vpsq(feat_id);
                set_value(pbc_matrix, feat_id, cell_id,
                    value_at(pbc_matrix, feat_id, cell_id) * vpsq(feat_id));
            }
        }
        // log_mem_usage('Updated Batch Adjustments')
        // for i in range(len(stdcfg_counts)):
        //     bayesdata[:, stdcfg_lookup == i] += stand_mean[:, i].reshape(
        //         (n_features, 1))
        // _bcdata.brief_print();
#pragma omp parallel for
        for (auto cell_id = 0; cell_id < _ncells; cell_id++) {
            IntType ex = _ebt_config.lookup(cell_id);
#pragma omp parallel for
            for (auto feat_id = 0; feat_id < _nfeatures; feat_id++) {
                // (*pbc_matrix)(feat_id, cell_id) += _stand_mean(feat_id, ex);
                set_value(pbc_matrix, feat_id, cell_id,
                    value_at(pbc_matrix, feat_id, cell_id) +
                    _stand_mean(feat_id, ex));
            }
        }
        // log_mem_usage('Normalized data')
        return pbc_matrix;
    }

  public:
    explicit ScementImpl(const SpMatType& spx, const FMatType& dsgn,
                         const FMatType& bdsgn, const SparseCfgType& cfg,
                         const SparseCfgType& ecfg, const FMatType& bcts,
                         IndexType nctx, IndexType nbtx)
        : _nthreads(omp_get_max_threads()), _adata(spx), _design(dsgn),
          _batch_design(bdsgn), _config(cfg), _ebt_config(ecfg),
          _batch_counts(bcts), _ncovariates(nctx), _nbatches(nbtx),
          _nconfig(cfg.n_config), _d_prime(_adata, _adata.data),
          _nfeatures(_adata.n_rows), _ncells(_adata.n_cols),
          _nebtcfg(_ebt_config.n_config), _cell_select(_ncells),
          _feat_select(_nfeatures), _cell_select_all(_ncells),
          _feat_select_all(_nfeatures), _delta_hat(_nbatches),
          _aprior(_nbatches), _bprior(_nbatches), _param_star(_nbatches),
          _batch_adj(_nbatches), _ebt_counts(_nebtcfg, _nbatches) {
        assert(_nbatches <= _ncovariates);
        assert(_nebtcfg <= _ncovariates);
        init();
    }

    template <typename PrintHelper>
    void adjust_be(const PrintHelper& prtr, std::ostream& ox, bool prt_debug,
                   bool prt_timings) {
        //
        timer run_timer;
        run_timer.reset();
        bhat();
        PRINT_IF(prt_debug, prtr.print(B_hat(), "B hat:", ox));
        PRINT_IF(prt_debug, prtr.print(cfgB_hat(), "Cfg B hat:", ox));
        PRINT_IF(prt_timings,
                 ox << "BHAT TIME (ms) :" << run_timer.elapsed() << std::endl);
        //
        run_timer.reset();
        compute_var_pooled();
        PRINT_IF(prt_timings, ox << "VAR POOLED TIME (ms) :"
                                   << run_timer.elapsed() << std::endl);
        PRINT_IF(prt_debug, prtr.print(var_pooled(), "var_pooled", ox));
        //
        stand_mean();
        PRINT_IF(prt_timings, ox << "STAND MEAN TIME (ms) :"
                                   << run_timer.elapsed() << std::endl);
        PRINT_IF(prt_debug, prtr.print(std_mean(), "stand_mean", ox));
        //
        run_timer.reset();
        standardize_data();
        PRINT_IF(prt_timings, ox << "STAND DAT TIME (ms) :"
                                   << run_timer.elapsed() << std::endl);
        PRINT_IF(prt_debug, prtr.print(m_prime(), "_m_prime", ox));
        PRINT_IF(prt_debug, prtr.print(d_prime(), "_d_prime", ox));
        //
        run_timer.reset();
        gamma_hat_prime();
        PRINT_IF(prt_timings, ox << "GAMAA HAT MEAN TIME (ms) :"
                                   << run_timer.elapsed() << std::endl);
        PRINT_IF(prt_debug, prtr.print(ghat_pfx(), "_ghat_pfx", ox));
        PRINT_IF(prt_debug, prtr.print(ghat_dprime(), "_ghat_dprime", ox));
        PRINT_IF(prt_debug, prtr.print(ghat_mprime(), "_ghat_mprime", ox));
        PRINT_IF(prt_debug, prtr.print(ghat(), "_ghat", ox));
        //
        run_timer.reset();
        multplicative_effect();
        PRINT_IF(prt_timings, ox << "MULT EFFECT TIME (ms) :"
                                   << run_timer.elapsed() << std::endl);
        PRINT_IF(
            prt_debug, for (int i = 0; i < delta_hat().size(); i++) {
                std::string batch_pfx = "_delta_hat " + std::to_string(i);
                prtr.print(delta_hat()[i], batch_pfx, ox);
            });
        //
        run_timer.reset();
        fix_hyperparams();
        PRINT_IF(prt_timings, ox << "HPARAM TIME (ms) :"
                                   << run_timer.elapsed() << std::endl);
        PRINT_IF(prt_debug, prtr.print(gamma_bar(), "_gamma_bar", ox));
        PRINT_IF(prt_debug, prtr.print(t2(), "_t2", ox));
        PRINT_IF(prt_debug, prtr.print(aprior(), "_aprior", ox));
        PRINT_IF(prt_debug, prtr.print(bprior(), "_bprior", ox));
        //
        run_timer.reset();
        parametric_adjustments();
        PRINT_IF(prt_timings, ox << "PARAM ADJ TIME (ms) :"
                                   << run_timer.elapsed() << std::endl);
        PRINT_IF(
            prt_debug, for (int i = 0; i < param_star().size(); i++) {
                std::string batch_pfx = "_param_star " + std::to_string(i);

                prtr.print(param_star()[i], batch_pfx, ox);
            });
        //
        run_timer.reset();
        adjust_data();
        PRINT_IF(prt_timings, ox << "ADJ DATA TIME (ms) :"
                                   << run_timer.elapsed() << std::endl);
        PRINT_IF(
            prt_debug, for (int i = 0; i < batch_adj().size(); i++) {
                std::string batch_pfx = "_batch_adj " + std::to_string(i);
                prtr.print(batch_adj()[i], batch_pfx, ox);
            });
    }
};

//
template <typename SCTImpl, typename DLImpl>
void run_bec(DLImpl* ploader, std::ostream& ox, bool print_debug = false,
             bool print_timings = false) {
    timer run_timer;
    SCTImpl sctimx(ploader);
    PRINT_IF(print_debug, sctimx.print(ox));
    PRINT_IF(print_timings,
             ox << "LOAD TIME (ms) :" << run_timer.elapsed() << std::endl);
    sctimx.impl().adjust_be(sctimx.prn(), ox, print_debug, print_timings);
    //
    run_timer.reset();
    typename SCTImpl::BECTMatType bcdata(sctimx.impl().nfeatures(),
                                         sctimx.impl().ncells());
    sctimx.impl().construct_bcdata(&bcdata);
    PRINT_IF(print_timings,
             ox << "BCDATA TIME (ms) :" << run_timer.elapsed() << std::endl);
    PRINT_IF(print_debug, sctimx.prn().print(bcdata, "_bcdata", ox));
}

#endif  // !IMPL_BASE_HPP
