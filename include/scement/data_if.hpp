#ifndef DATA_INTERFACE_HPP
#define DATA_INTERFACE_HPP

template<typename SDT>
class DataInterface {
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

public: // Member functions
    virtual FMatType batch_counts() = 0;
    virtual FMatType design() = 0;
    virtual SparseCfgType config() = 0;
    virtual SparseCfgType ebt_config() = 0;
    virtual SpMatType adata() = 0;
    virtual ~DataInterface() {}
};


#endif // !DATA_LOADER_HPP
