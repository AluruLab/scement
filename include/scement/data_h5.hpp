#ifndef HDF5_EIGEN_LOADER_HPP
#define HDF5_EIGEN_LOADER_HPP

// Hivefive
#include <Eigen/Dense>
#include <highfive/highfive.hpp>
#include <string>
#include <vector>
//
#include <scement/coocsr.hpp>
#include <scement/data_if.hpp>
#include <scement/impl_eigen.hpp>

template <typename IndexType>
inline static IndexType HF_vector_size(HighFive::DataSet h5data) {
    return IndexType(h5data.getDimensions().at(0));
}

template <typename IndexType>
inline static IndexType HF_matrix_rows(HighFive::DataSet h5data) {
    return IndexType(h5data.getDimensions().at(0));
}

template <typename IndexType>
inline static IndexType HF_matrix_columns(HighFive::DataSet h5data) {
    return IndexType(h5data.getDimensions().at(1));
}

template <typename IndexType>
static void HF_matrix_size(HighFive::DataSet h5data, IndexType* nrows,
                           IndexType* ncols) {
    auto dims = h5data.getDimensions();
    *nrows = dims[0];
    *ncols = dims[1];
}

template <typename PElemType>
static void HF_get_pair(HighFive::DataSet h5data, PElemType* px,
                        PElemType* py) {
    std::vector<PElemType> vpairs;
    h5data.read(vpairs);
    *px = vpairs[0];
    *py = vpairs[1];
}


struct HF2Eigen {
    template <typename VecType, typename VElemType>
    inline static VecType vec(HighFive::DataSet h5data) {
        auto dims = h5data.getDimensions();
        VecType result(dims[0]);
        h5data.read<VElemType>(result.data());
        return result;
    }

    template <typename MatType, typename MElemType>
    inline static MatType mat_transposed(HighFive::DataSet h5data) {
        auto dims = h5data.getDimensions();
        MatType result(dims[1], dims[0]);
        h5data.read<MElemType>(result.data());
        return result;
    }

    template <typename MatType, typename MElemType>
    inline static MatType mat(HighFive::DataSet h5data) {
        auto dims = h5data.getDimensions();
        MatType result(dims[1], dims[0]);
        h5data.read<MElemType>(result.data());
        return result.transpose();
    }

    template <typename MatType, typename VElemType>
    inline static MatType vec2mat(HighFive::DataSet h5data) {
        auto dims = h5data.getDimensions();
        MatType result(dims[0], 1);
        h5data.read<VElemType>(result.data());
        return result;
    }
};

template<typename SDT>
class H5EigenDIaT : public DataInterface<SDT> {
public:
    using ElemType = typename SDT::ElemType;
    using IntType = typename SDT::IntType;
    using IndexType = typename SDT::IndexType;
    //
    using FMatType = typename SDT::FMatType;
    using IMatType = typename SDT::IMatType;
    using IVecType = typename SDT::IVecType;
    //
    using SparseCfgType = typename SDT::SparseCfgType;
    using SpMatType = typename SDT::SpMatType;

private:
    const HighFive::File& _h5fx;
    IndexType _ngenes, _ncells, _nelems;
    IndexType _nconfig, _nebtcfg;

  public:  // Member Functions
    inline SpMatType adata() { 
        return SpMatType(
            _ngenes, _ncells, _nelems,
            HF2Eigen::vec2mat<FMatType, ElemType>(_h5fx.getDataSet("/uns/data")),
            HF2Eigen::mat<IMatType, IntType>(_h5fx.getDataSet("/uns/coo")),
            HF2Eigen::mat<IMatType, IntType>(_h5fx.getDataSet("/uns/csrptr")));
    }

    inline SparseCfgType config() {
        return SparseCfgType(
            HF2Eigen::mat<IMatType, int32_t>(_h5fx.getDataSet("/uns/cuniq")),
            HF2Eigen::vec<IVecType, int32_t>(_h5fx.getDataSet("/uns/ccounts")),
            HF2Eigen::vec<IVecType, int32_t>(_h5fx.getDataSet("/uns/clookup")),
            _nconfig);
    }

    inline SparseCfgType ebt_config() {
        return SparseCfgType(
            HF2Eigen::mat<IMatType, int32_t>(_h5fx.getDataSet("/uns/euniq")),
            HF2Eigen::vec<IVecType, int32_t>(_h5fx.getDataSet("/uns/ecounts")),
            HF2Eigen::vec<IVecType, int32_t>(_h5fx.getDataSet("/uns/elookup")),
            _nebtcfg);
    }

    inline FMatType design() {
        return HF2Eigen::mat<FMatType, float>(_h5fx.getDataSet("/uns/design"));
    }

    inline FMatType batch_counts() {
        return HF2Eigen::mat<FMatType, float>(
            _h5fx.getDataSet("/uns/batch_counts"));
    }

    explicit H5EigenDIaT(const HighFive::File& h5file) : _h5fx(h5file) {
        _nconfig = HF_matrix_rows<IndexType>(_h5fx.getDataSet("/uns/cuniq"));
        _nebtcfg = HF_matrix_rows<IndexType>(_h5fx.getDataSet("/uns/euniq"));
        _nelems = HF_vector_size<IndexType>(_h5fx.getDataSet("/uns/data"));
        HF_get_pair<IndexType>(_h5fx.getDataSet("/uns/shape"), &_ngenes,
                               &_ncells);
    }
};

using H5EigenDI = H5EigenDIaT<ScementEigen<>>;

#ifdef SCEMENT_ARMA
// Armadillo
#include <scement/impl_arma.hpp>
#include <armadillo>

struct HF2Arma {
    template <typename VecType, typename VElemType>
    inline static VecType vec(HighFive::DataSet ldataset) {
        auto dims = ldataset.getDimensions();
        VecType result(dims[0]);
        ldataset.read<VElemType>(result.memptr());
        return result;
    }

    template <typename MatType, typename VElemType>
    inline static MatType vec2mat(HighFive::DataSet h5data) {
        auto dims = h5data.getDimensions();
        MatType result(dims[0], 1);
        h5data.read<VElemType>(result.memptr());
        return result;
    }

    template <typename MatType, typename MElemType>
    inline static MatType mat_transposed(HighFive::DataSet h5data) {
        auto dims = h5data.getDimensions();
        MatType result(dims[1], dims[0]);
        h5data.read<MElemType>(result.memptr());
        return result;
    }

    template <typename MatType, typename MElemType>
    inline static MatType mat(HighFive::DataSet h5data) {
        auto dims = h5data.getDimensions();
        MatType result(dims[1], dims[0]);
        h5data.read<MElemType>(result.memptr());
        return result.t();
    }
};

template<typename SDT>
class H5ArmaDIaT : public DataInterface<SDT> {
public: // Types
    using ElemType = typename SDT::ElemType;
    using IntType = typename SDT::IntType;
    using IndexType = typename SDT::IndexType;
    //
    using FMatType = typename SDT::FMatType;
    using IMatType = typename SDT::IMatType;
    using IVecType = typename SDT::IVecType;
    //
    using SparseCfgType = typename SDT::SparseCfgType;
    using SpMatType = typename SDT::SpMatType;
private:
    const HighFive::File& _h5fx;
    IndexType _ngenes, _ncells, _nelems, _nconfig, _nebtcfg;
public: // Member Functions
    virtual inline SpMatType adata() { 
        return SpMatType(
            _ngenes, _ncells, _nelems,
            HF2Arma::vec2mat<FMatType, ElemType>(_h5fx.getDataSet("/uns/data")),
            HF2Arma::mat<IMatType, IntType>(_h5fx.getDataSet("/uns/coo")),
            HF2Arma::mat<IMatType, IntType>(_h5fx.getDataSet("/uns/csrptr")));
    }

    virtual inline SparseCfgType config() {
        return SparseCfgType(
            HF2Arma::mat<IMatType, int32_t>(_h5fx.getDataSet("/uns/cuniq")),
            HF2Arma::vec<IVecType, int32_t>(_h5fx.getDataSet("/uns/ccounts")),
            HF2Arma::vec<IVecType, int32_t>(_h5fx.getDataSet("/uns/clookup")),
            _nconfig);
    }

    virtual inline SparseCfgType ebt_config() {
        return SparseCfgType(
            HF2Arma::mat<IMatType, int32_t>(_h5fx.getDataSet("/uns/euniq")),
            HF2Arma::vec<IVecType, int32_t>(_h5fx.getDataSet("/uns/ecounts")),
            HF2Arma::vec<IVecType, int32_t>(_h5fx.getDataSet("/uns/elookup")),
            _nebtcfg);
    }

    virtual inline FMatType design() {
        return HF2Arma::mat<FMatType, float>(_h5fx.getDataSet("/uns/design"));
    }

    virtual inline FMatType batch_counts() {
        return HF2Arma::mat<FMatType, float>(
            _h5fx.getDataSet("/uns/batch_counts"));
    }

    explicit H5ArmaDIaT(const HighFive::File& h5file) : _h5fx(h5file) {
        _nconfig = HF_matrix_rows<IndexType>(_h5fx.getDataSet("/uns/cuniq"));
        _nebtcfg = HF_matrix_rows<IndexType>(_h5fx.getDataSet("/uns/euniq"));
        _nelems = HF_vector_size<IndexType>(_h5fx.getDataSet("/uns/data"));
        HF_get_pair<IndexType>(_h5fx.getDataSet("/uns/shape"), &_ngenes,
                               &_ncells);
    }
};

using H5ArmaDI = H5ArmaDIaT<ScementArma<>>;

#endif

#endif  // !HDF5_EIGEN_LOADER_HPP
