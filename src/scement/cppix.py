from typing import Literal, List
import pandas as pd
import pathlib
import anndata
import numpy as np

#
from scipy.sparse import csr_matrix
from scanpy import logging, settings
from . import log_mem_usage, sparse_design_matrix
from ._scementcpp import bec_eg, bec_arma, arma_available


class ScementCPP:
    def __init__(
        self, axcat: anndata.AnnData, adata: csr_matrix, key: str, cvats: list[str]
    ):
        self.axcat: anndata.AnnData = axcat
        self.adata: csr_matrix = adata
        self.batch_key = key
        self.covariates = cvats
        self.astyle: Literal["K", "A", "C", "F"] = "F"
        self.ftype = np.float32
        self.itype = np.int32
        self.rdict = {}
        logging._set_log_level(settings, logging.INFO)

    def design_config(self):
        # construct a pandas series of the batch annotation
        obs_colums = [self.batch_key] + list(self.covariates)
        self.model: pd.DataFrame = self.axcat.obs[obs_colums]  # type:ignore
        self.batch_info = self.model.groupby(self.batch_key).indices.values()
        batch_items = self.model.groupby(self.batch_key, observed=False).groups.items()
        self.batch_levels, self.batch_info = zip(*batch_items)
        self.n_batch = len(self.batch_info)
        self.n_batches = np.array([len(v) for v in self.batch_info])
        self.n_features = self.adata.shape[0]
        log_mem_usage("Model Construted")
        # Design Matrix
        self.design = sparse_design_matrix(
            self.model,  # type:ignore
            self.batch_key,
            self.batch_levels,
        )
        log_mem_usage("Design Completed: " + str(self.design.shape))
        # Design Config : covariates
        self.icfg_vectors, self.cfg_lookup, self.cfg_counts = np.unique(
            self.design, axis=0, return_counts=True, return_inverse=True
        )
        # Design Config : extra-batch variables
        self.tmp: np.ndarray = np.array(self.design.copy())
        self.tmp[:, : self.n_batch] = 0
        self.istdcfg_vecs, self.stdcfg_lookup, self.stdcfg_counts = np.unique(
            self.tmp, axis=0, return_counts=True, return_inverse=True
        )
        log_mem_usage("Formula Completed: " + str(self.design.shape))
        #
        adata_coo = self.adata.tocoo()
        # CSRPTR and COO
        rptr_shape = (len(self.adata.indptr) - 1, 2)
        self.csrptr = np.zeros(
            rptr_shape, dtype=self.itype, order=self.astyle
        )  # type:ignore
        self.csrptr[:, 0] = self.adata.indptr[:-1]
        self.csrptr[:, 1] = self.adata.indptr[1:]
        coo_shape = (len(adata_coo.row), 2)
        self.coo = np.zeros(
            coo_shape, dtype=self.itype, order=self.astyle
        )  # type:ignore
        self.coo[:, 0] = adata_coo.row
        self.coo[:, 1] = adata_coo.col
        # Batch Counts
        self.batch_counts = np.array(
            # np.array(self.n_batches, ndmin=2, dtype=self.ftype).T,
            self.n_batches,
            order=self.astyle,
            dtype=self.ftype,
        ).reshape((len(self.n_batches), 1))
        log_mem_usage("PREP Completed: " + str(self.design.shape))
        return self

    def prep_dict(self):
        self.rdict["ngenes"] = np.uint32(self.adata.shape[0])
        self.rdict["ncells"] = np.uint32(self.adata.shape[1])
        self.rdict["nelems"] = np.uint32(self.adata.data.shape[0])
        self.rdict["shape"] = np.array(
            self.adata.shape, dtype=self.itype, order=self.astyle
        )
        # COOCSR data array
        # in_dict["data"] = np.array(self.adata.data, dtype=np.float32)
        self.rdict["data"] = self.adata.data
        self.rdict["coo"] = self.coo
        self.rdict["csrptr"] = self.csrptr
        self.rdict["batch_counts"] = self.batch_counts
        # Design
        # self.np_design = np.array(self.design, dtype=ftype, order=self.astyle)
        # self.rdict["design"] = self.np_design
        self.rdict["design"] = np.array(
            self.design.to_numpy(copy=False), dtype=self.ftype, order=self.astyle
        )
        # Design covariates config
        self.rdict["cuniq"] = np.array(
            self.icfg_vectors, dtype=self.itype, order=self.astyle
        )
        self.rdict["clookup"] = np.array(
            self.cfg_lookup, dtype=self.itype, order=self.astyle
        )
        self.rdict["ccounts"] = np.array(
            self.cfg_counts.reshape(self.cfg_counts.size, 1),
            dtype=self.itype,
            order=self.astyle,
        )
        # Design extra-batch config
        self.rdict["euniq"] = np.array(
            self.istdcfg_vecs, dtype=self.itype, order=self.astyle
        )
        self.rdict["elookup"] = np.array(
            self.stdcfg_lookup, dtype=self.itype, order=self.astyle
        )
        self.rdict["ecounts"] = np.array(
            self.stdcfg_counts.reshape(self.stdcfg_counts.size, 1),
            dtype=self.itype,
            order=self.astyle,
        )
        return self

    def prep(self, out_file: str | None = None):
        if not self.rdict:
            self.design_config()
            self.prep_dict()
        if out_file is not None:
            for kx, vx in self.rdict.items():
                self.axcat.uns[kx] = vx
            self.axcat.write_h5ad(pathlib.PurePath(out_file))
        return self

    def run(
        self,
        axcat: anndata.AnnData,
        inplace: bool,
        print_debug: bool,
        print_timings: bool,
    ):
        if not self.rdict:
            self.design_config()
            self.prep_dict()
        if arma_available():
            bayesdata = bec_arma(self.rdict, print_debug, print_timings)
        else:
            bayesdata = bec_eg(self.rdict, print_debug, print_timings)
        if inplace:
            axcat.X = bayesdata
        else:
            return bayesdata

    @classmethod
    def from_anndata(cls, axcat: anndata.AnnData, key: str, cvats: List[str]):
        log_mem_usage("START")
        axcat._sanitize()
        adata = axcat.X.transpose().tocsr()  # type:ignore
        log_mem_usage("Transposed Input Matrix")
        return cls(axcat, adata, key, cvats)

    @classmethod
    def from_file(cls, input_file: str, key: str, cvats: List[str]):
        log_mem_usage("START")
        if input_file.endswith("loom"):
            axcat = anndata.read_loom(pathlib.PurePath(input_file))
        else:
            axcat = anndata.read_h5ad(input_file)
        axcat._sanitize()
        log_mem_usage("Loaded " + input_file)
        adata = axcat.X.transpose().tocsr()  # type:ignore
        log_mem_usage("Transposed Input Matrix")
        return cls(axcat, adata, key, cvats)


def prep_anndata_cpp(
    input_file: str, key: str, cvats: List[str], output_file: str | None = None
):
    spx = ScementCPP.from_file(input_file, key, cvats)
    spx.prep(output_file)
    return spx


def sct_sparse_cpp(
    axcat: anndata.AnnData,
    key: str,
    cvats: List[str],
    inplace: bool = True,
    print_debug: bool = False,
    print_timings: bool = False,
):
    spx = ScementCPP.from_anndata(axcat, key, cvats)
    spx.run(axcat, inplace, print_debug, print_timings)
