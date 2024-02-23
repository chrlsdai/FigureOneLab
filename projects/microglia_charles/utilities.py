import scipy
import scanpy as sc
import anndata


def sampleQC(
    adata, min_genes=200, percent_cells_threshold=0.1, percent_mt_threshold=10
):
    n_cells, n_genes = adata.shape
    print("Before QC: {} cells, {} genes".format(n_cells, n_genes))

    sc.pp.filter_cells(adata, min_genes=min_genes)
    if percent_cells_threshold:
        sc.pp.filter_genes(
            adata, min_cells=round(percent_cells_threshold * n_cells / 100)
        )

    # adata.layers["counts"] = scipy.sparse.csr_matrix(adata.X.copy())
    # adata.layers['log2_counts'] = scipy.sparse.csr_matrix(sc.pp.log1p(adata.layers['counts'],copy(), base=2))

    if percent_mt_threshold:
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )
        adata.obs.rename(columns={"pct_counts_mt": "percent_mt"}, inplace=True)
        adata.var.drop("mt", axis=1, inplace=True)
        adata = adata[adata.obs["percent_mt"] < percent_mt_threshold,].copy()

    n_cells, n_genes = adata.shape
    print("After QC: {} cells, {} genes".format(n_cells, n_genes))

    return adata
