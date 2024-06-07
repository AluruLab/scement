//
// Copyright [2024]
//
#ifndef PY_LOADER_HPP
#define PY_LOADER_HPP

#include <Eigen/Core>
#include <scement/impl_eigen.hpp>
#include <scement/py_types.hpp>

template <typename IndexType>
static inline IndexType py_config_size(pydict in_dict, const char* uniq_key) {
    pyiarray_t py_unq = in_dict[uniq_key].cast<pyiarray_t>();
    return IndexType(py_unq.shape(0));
}

//
struct Numpy2Eigen {

    template <typename PyType, typename MapMatType, typename MatElemType>
    static inline MapMatType mapmat(pydict in_dict, const char* dkey) {
        PyType pymat = in_dict[dkey].cast<PyType>();
        // MatElemType* fptr_mat =
        // static_cast<MatElemType*>(pymat.request().ptr);
        return MapMatType(static_cast<MatElemType*>(pymat.request().ptr),
                          pymat.shape(0), pymat.shape(1));
    }

    template <typename PyType, typename MatType, typename MatElemType>
    static inline MatType mat(pydict in_dict, const char* dkey) {
        return mapmat<PyType, Eigen::Map<MatType>, MatElemType>(in_dict, dkey);
    }

    template <typename PyType, typename MapMatType, typename MatElemType>
    static inline MapMatType vec2mapmat(pydict in_dict, const char* dkey) {
        PyType pymat = in_dict[dkey].cast<PyType>();
        // std::size_t nrows = pymat.shape(0);
        // pybuffer_info buf_mat = pymat.request();
        // MatElemType* fptr_mat = static_cast<MatElemType*>(buf_mat.ptr);
        return MapMatType(static_cast<MatElemType*>(pymat.request().ptr),
                          pymat.shape(0), 1);
    }

    template <typename PyType, typename MatType, typename MatElemType>
    static inline MatType vec2mat(pydict in_dict, const char* dkey) {
        return vec2mapmat<PyType, Eigen::Map<MatType>, MatElemType>(in_dict,
                                                                    dkey);
    }

    template <typename PyType, typename MapVecType, typename MatElemType>
    static inline MapVecType vec2mapvec(pydict in_dict, const char* dkey) {
        PyType pymat = in_dict[dkey].cast<PyType>();
        // std::size_t nrows = pymat.shape(0);
        // pybuffer_info buf_mat = pymat.request();
        // MatElemType* fptr_mat = static_cast<MatElemType*>(buf_mat.ptr);
        //
        return MapVecType(static_cast<MatElemType*>(pymat.request().ptr),
                          pymat.shape(0));
    }

    template <typename PyType, typename VecType, typename MatElemType>
    static inline VecType vec2vec(pydict in_dict, const char* dkey) {
        return vec2mapmat<PyType, Eigen::Map<VecType>, MatElemType>(in_dict,
                                                                    dkey);
    }

};

template <typename SDT> class NumPyEigenDIaT : public DataInterface<SDT> {
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
    const pydict& _data_dict;
    IndexType _nfeatures, _ncells, _nelems, _nconfig, _nebtcfg;

  public:  // Member Functions
    //
    virtual inline SpMatType adata() {
        std::cout << "NGENES   : " << _nfeatures << std::endl
                  << "NCELLS   : " << _ncells << std::endl
                  << "NCOVAR   : " << _nelems << std::endl;
        return SpMatType(
            _nfeatures, _ncells, _nelems,
            Numpy2Eigen::vec2mat<pyfarray_t, FMatType, float>(_data_dict,
                                                              "data"),
            Numpy2Eigen::mat<pyiarray_t, IMatType, int32_t>(_data_dict, "coo"),
            Numpy2Eigen::mat<pyiarray_t, IMatType, int32_t>(_data_dict,
                                                            "csrptr"));
    }

    virtual inline SparseCfgType config() {
        return SparseCfgType(
            Numpy2Eigen::mat<pyiarray_t, IMatType, int32_t>(_data_dict,
                                                            "cuniq"),
            Numpy2Eigen::mat<pyiarray_t, IMatType, int32_t>(_data_dict,
                                                            "ccounts"),
            Numpy2Eigen::vec2vec<pyiarray_t, IVecType, int32_t>(_data_dict,
                                                                "clookup"),
            _nconfig);
    }

    virtual inline SparseCfgType ebt_config() {
        return SparseCfgType(
            Numpy2Eigen::mat<pyiarray_t, IMatType, int32_t>(_data_dict,
                                                            "euniq"),
            Numpy2Eigen::mat<pyiarray_t, IMatType, int32_t>(_data_dict,
                                                            "ecounts"),
            Numpy2Eigen::vec2vec<pyiarray_t, IVecType, int32_t>(_data_dict,
                                                                "elookup"),
            _nebtcfg);
    }

    virtual inline FMatType design() {
        return Numpy2Eigen::mat<pyfarray_t, FMatType, float>(_data_dict,
                                                             "design");
    }

    virtual inline FMatType batch_counts() {
        return Numpy2Eigen::mat<pyfarray_t, FMatType, float>(_data_dict,
                                                             "batch_counts");
    }

    explicit NumPyEigenDIaT(const pydict& pdict)
        : _data_dict(pdict), _nfeatures(_data_dict["ngenes"].cast<IndexType>()),
          _ncells(_data_dict["ncells"].cast<IndexType>()),
          _nelems(_data_dict["nelems"].cast<IndexType>()),
          _nconfig(py_config_size<IndexType>(_data_dict, "cuniq")),
          _nebtcfg(py_config_size<IndexType>(_data_dict, "euniq")) {}
};

template <typename SDT> class NumPyEigenMapDIaT : public DataInterface<SDT> {
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
    //
    using MapIMatType = typename SDT::MapIMatType;
    using MapFMatType = typename SDT::MapFMatType;
    using MapIVecType = typename SDT::MapIVecType;
    //
    const pydict& _data_dict;
    IndexType _nfeatures, _ncells, _nelems, _nconfig, _nebtcfg;

  public:
    //
    virtual inline SpMatType adata() {
        // std::cout << "NGENES   : " << _nfeatures << std::endl
        //           << "NCELLS   : " << _ncells << std::endl
        //           << "NCOVAR   : " << _nelems << std::endl;
        return SpMatType(
            _nfeatures, _ncells, _nelems,
            Numpy2Eigen::vec2mapmat<pyfarray_t, MapFMatType, ElemType>(
                _data_dict, "data"),
            Numpy2Eigen::mapmat<pyiarray_t, MapIMatType, IntType>(_data_dict,
                                                                  "coo"),
            Numpy2Eigen::mapmat<pyiarray_t, MapIMatType, IntType>(_data_dict,
                                                                  "csrptr"));
    }

    virtual inline SparseCfgType config() {
        return SparseCfgType(
            Numpy2Eigen::mapmat<pyiarray_t, MapIMatType, IntType>(_data_dict,
                                                                  "cuniq"),
            Numpy2Eigen::mapmat<pyiarray_t, MapIMatType, IntType>(_data_dict,
                                                                  "ccounts"),
            Numpy2Eigen::vec2mapvec<pyiarray_t, MapIVecType, IntType>(
                _data_dict, "clookup"),
            _nconfig);
    }

    virtual inline SparseCfgType ebt_config() {
        return SparseCfgType(
            Numpy2Eigen::mapmat<pyiarray_t, MapIMatType, IntType>(_data_dict,
                                                                  "euniq"),
            Numpy2Eigen::mapmat<pyiarray_t, MapIMatType, IntType>(_data_dict,
                                                                  "ecounts"),
            Numpy2Eigen::vec2mapvec<pyiarray_t, MapIVecType, IntType>(
                _data_dict, "elookup"),
            _nebtcfg);
    }

    virtual inline FMatType design() {
        return Numpy2Eigen::mapmat<pyfarray_t, MapFMatType, ElemType>(
            _data_dict, "design");
    }

    virtual inline FMatType batch_counts() {
        return Numpy2Eigen::mapmat<pyfarray_t, MapFMatType, ElemType>(
            _data_dict, "batch_counts");
    }

    explicit NumPyEigenMapDIaT(const pydict& pdict)
        : _data_dict(pdict), _nfeatures(_data_dict["ngenes"].cast<IndexType>()),
          _ncells(_data_dict["ncells"].cast<IndexType>()),
          _nelems(_data_dict["nelems"].cast<IndexType>()),
          _nconfig(py_config_size<IndexType>(_data_dict, "cuniq")),
          _nebtcfg(py_config_size<IndexType>(_data_dict, "euniq")) {}
};

#ifdef SCEMENT_ARMA
#include <scement/impl_arma.hpp>

struct Numpy2Arma {

    template <typename PyType, typename MatType, typename ElemType>
    static inline MatType mat(pydict in_dict, const char* dkey) {
        PyType pymat = in_dict[dkey].cast<PyType>();
        // ElemType* fptr_mat = static_cast<MatElemType*>(pymat.request().ptr);
        // return MatType(static_cast<ElemType*>(pymat.request().ptr),
        //               pymat.shape(0), pymat.shape(1), copy_aux_mem = false);
        return MatType(static_cast<ElemType*>(pymat.request().ptr),
                       pymat.shape(0), pymat.shape(1), false, false);
    }

    template <typename PyType, typename VecType, typename ElemType>
    static inline VecType vec(pydict in_dict, const char* dkey) {
        PyType pymat = in_dict[dkey].cast<PyType>();
        bool copy_aux_mem = false;
        // ElemType* fptr_mat = static_cast<MatElemType*>(pymat.request().ptr);
        // return VecType(static_cast<ElemType*>(pymat.request().ptr),
        //              pymat.shape(0), copy_aux_mem = false);
        return VecType(static_cast<ElemType*>(pymat.request().ptr),
                       pymat.shape(0), false, false);
    }

    template <typename PyType, typename MatType, typename ElemType>
    static inline MatType vec2mat(pydict in_dict, const char* dkey) {
        PyType pymat = in_dict[dkey].cast<PyType>();
        // bool copy_aux_mem;
        //  pybuffer_info buf_mat = pymat.request();
        //  ElemType* fptr_mat = static_cast<ElemType*>(buf_mat.ptr);
        //  return MatType(static_cast<ElemType*>(pymat.request().ptr),
        //                pymat.shape(0), 1, copy_aux_mem = false);
        return MatType(static_cast<ElemType*>(pymat.request().ptr),
                       pymat.shape(0), 1, false, false);
    }
};

template <typename SDT> class NumPyArmaDIaT : public DataInterface<SDT> {
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
    const pydict& _data_dict;
    IndexType _nfeatures, _ncells, _nelems;
    IndexType _nconfig, _nebtcfg;

  public:  // Member Functions
    //
    virtual inline SparseCfgType config() {
        return SparseCfgType(
            Numpy2Arma::mat<pyiarray_t, IMatType, IntType>(_data_dict, "cuniq"),
            Numpy2Arma::mat<pyiarray_t, IMatType, IntType>(_data_dict,
                                                           "ccounts"),
            Numpy2Arma::vec<pyiarray_t, IVecType, IntType>(_data_dict,
                                                           "clookup"),
            _nconfig);
    }

    virtual inline SpMatType adata() {
        // std::cout << "NGENES   : " << _nfeatures << std::endl
        //           << "NCELLS   : " << _ncells << std::endl
        //           << "NCOVAR   : " << _nelems << std::endl;
        return SpMatType(
            _nfeatures, _ncells, _nelems,
            Numpy2Arma::vec2mat<pyfarray_t, FMatType, ElemType>(_data_dict,
                                                                "data"),
            Numpy2Arma::mat<pyiarray_t, IMatType, IntType>(_data_dict, "coo"),
            Numpy2Arma::mat<pyiarray_t, IMatType, IntType>(_data_dict,
                                                           "csrptr"));
    }

    virtual inline SparseCfgType ebt_config() {
        return SparseCfgType(
            Numpy2Arma::mat<pyiarray_t, IMatType, IntType>(_data_dict, "euniq"),
            Numpy2Arma::mat<pyiarray_t, IMatType, IntType>(_data_dict,
                                                           "ecounts"),
            Numpy2Arma::vec<pyiarray_t, IVecType, IntType>(_data_dict,
                                                           "elookup"),
            _nebtcfg);
    }

    virtual inline FMatType design() {
        return Numpy2Arma::mat<pyfarray_t, FMatType, ElemType>(_data_dict,
                                                               "design");
    }

    virtual inline FMatType batch_counts() {
        return Numpy2Arma::mat<pyfarray_t, FMatType, ElemType>(_data_dict,
                                                               "batch_counts");
    }

    explicit NumPyArmaDIaT(const pydict& pdict)
        : _data_dict(pdict), _nfeatures(_data_dict["ngenes"].cast<IndexType>()),
          _ncells(_data_dict["ncells"].cast<IndexType>()),
          _nelems(_data_dict["nelems"].cast<IndexType>()),
          _nconfig(py_config_size<IndexType>(_data_dict, "cuniq")),
          _nebtcfg(py_config_size<IndexType>(_data_dict, "euniq")) {}
};

#endif

#endif  //! PY_LOADER_HPP
