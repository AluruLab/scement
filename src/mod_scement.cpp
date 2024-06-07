//
// Copright [2024] <Sriram P C>
//
#include <scement/py_types.hpp>
#include <scement/impl_arma.hpp>
#include <scement/impl_eigen.hpp>
#include <scement/data_py.hpp>

using NumPyMatRef = pybind11::detail::unchecked_mutable_reference<float, 2>;
using NumPyEigenMapDI = NumPyEigenMapDIaT<ScementEigenMap<NumPyMatRef, true>>;
using NumPyArmaDI = NumPyArmaDIaT<ScementArma<NumPyMatRef, true>>;

void print_dict(pydict in_dict) {
    const char* dct_keys[]{"data",    "coo",     "csrptr",  "batch_counts",
                           "design",  "ccounts", "clookup", "cuniq",
                           "ecounts", "elookup", "euniq"};
    const char* float_array_keys[]{"data", "batch_counts", "design"};
    const char* int_array_keys[]{"coo",   "csrptr",  "ccounts", "clookup",
                                 "cuniq", "ecounts", "elookup", "euniq"};
    for (const char* pkey : float_array_keys) {
        pyfarray_t px = in_dict[pkey].cast<pyfarray_t>();
        pybind11::print(pkey, px.ndim(), px.shape(0),
                        px.ndim() > 1 ? px.shape(1) : 0);
    }
    for (const char* pkey : int_array_keys) {
        pyiarray_t px = in_dict[pkey].cast<pyiarray_t>();
        pybind11::print(pkey, px.ndim(), px.size(), px.shape(0),
                        px.ndim() > 1 ? px.shape(1) : 0);
    }
}

pyfarray_t bec_eg(pydict in_dict, bool print_debug = false,
                  bool print_timings = false) {
    timer run_timer;
    //
    NumPyEigenMapDI peg_ldr(in_dict);
    ScementEigenMap<NumPyMatRef, true> scteg(&peg_ldr);
    PRINT_IF(print_debug, scteg.print());
    //
    scteg.impl().adjust_be(scteg.prn(), std::cout, print_debug, print_timings);
    //
    // pyfarray_t bc_matrix({scteg.impl().nfeatures(), scteg.impl().ncells()});
    pyfarray_t bc_matrix({scteg.impl().ncells(), scteg.impl().nfeatures()});
    auto result = bc_matrix.mutable_unchecked<2>();
    scteg.impl().construct_bcdata(&result);
    PRINT_IF(print_timings,
             std::cout << "TOTAL TIME (ms) :" << run_timer.elapsed() << std::endl);
    // pyfarray_t bc_matrix({5, 12});
    return bc_matrix;
}

#ifdef SCEMENT_ARMA
pyfarray_t bec_arma(pydict in_dict, bool print_debug = false,
                    bool print_timings = false) {
    timer run_timer;
    NumPyArmaDI pa_loader(in_dict);
    ScementArma<NumPyMatRef, true> parma(&pa_loader);
    PRINT_IF(print_debug, parma.print());
    //
    parma.impl().adjust_be(parma.prn(), std::cout, print_debug, print_timings);
    //
    // pyfarray_t bc_matrix({parma.impl().nfeatures(), parma.impl().ncells()});
    pyfarray_t bc_matrix({parma.impl().ncells(), parma.impl().nfeatures()});
    auto result = bc_matrix.mutable_unchecked<2>();
    parma.impl().construct_bcdata(&result);
    PRINT_IF(print_timings,
             std::cout << "TOTAL TIME (ms) :" << run_timer.elapsed() << std::endl);
    return bc_matrix;
}

bool arma_available() { return true; }
#else
pyfarray_t bec_armax(pydict in_dict, bool transposed = false,
                    bool print_debug = false, bool print_timings = false) {
    std::cout << "Armadillo library NOT AVAILABLE" << std::endl;
    return pybind11::none();
}
bool arma_available() { return false; }
#endif

PYBIND11_MODULE(_scementcpp, m) {
    pyoptions options;
    // options.disable_function_signatures();

    m.doc() = "SCEMENT";
    m.def("bec_eg", &bec_eg, R"(
    Run SCEMENT Algorithm
    Args:
        in_dict dict with the following keys:
            data             (n_elem,)               np.ndarray[float32]
            coo              (n_elem, 2)             np.ndarray[int32]
            csrptr           (n_features, 2)         np.ndarray[int32] 
            ---
            batch_counts     (n_batches, 1)          np.ndarray[float32]
            design           (n_covar, n_features)   np.ndarray[float32]
            ---
            ccounts          (n_covar,)              np.ndarray[int32]
            clookup          (n_cells,)              np.ndarray[int32]
            cuniq            (n_covar, n_covar)      np.ndarray[int32]
            ---
            ecounts          (n_ebtcfg,)             np.ndarray[int32]
            elookup          (n_cells,)              np.ndarray[int32]
            euniq            (n_covar, n_ebtcfg)     np.ndarray[int32]
            ---
            shape            (2,)                    int32
            ncells           1                       uint32
            nfeatures        1                       uint32
            nelems           1                       uint32
          print_debug flag to print  debug statements default:False
          print_timings flag to print  debug statements default:False
    Returns:
        Batch Corrected Matrix )",
          pyarg("data_dict"), pyarg("print_debug"), pyarg("print_timings"));
    m.def("bec_arma", &bec_arma, R"(
    Run SCEMENT Algorithm
     Args:
        data_dict dict with the following keys:
            data             (n_elem,)               np.ndarray[float32]
            coo              (n_elem, 2)             np.ndarray[int32]
            csrptr           (n_features, 2)         np.ndarray[int32] 
            ---
            batch_counts     (n_batches, 1)          np.ndarray[float32]
            design           (n_covar, n_features)   np.ndarray[float32]
            ---
            ccounts          (n_covar,)              np.ndarray[int32]
            clookup          (n_cells,)              np.ndarray[int32]
            cuniq            (n_covar, n_covar)      np.ndarray[int32]
            ---
            ecounts          (n_ebtcfg,)             np.ndarray[int32]
            elookup          (n_cells,)              np.ndarray[int32]
            euniq            (n_covar, n_ebtcfg)     np.ndarray[int32]
            ---
            shape            (2,)                    int32
            ncells           1                       uint32
            nfeatures        1                       uint32
            nelems           1                       uint32
          print_debug flag to print  debug statements default:False
          print_timings flag to print  debug statements default:False
    Returns:
         Result)",
          pyarg("data_dict"), pyarg("print_debug"), pyarg("print_timings"));

    m.def("arma_available", &arma_available, "Is Armadillo available ?");
#ifdef SCEMENT_VERSION
    m.attr("__version__") = SCEMENT_VERSION;
#else
    m.attr("__version__") = "dev";
#endif
}
