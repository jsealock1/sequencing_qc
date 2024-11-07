# assign ancestry using gnomad v3 projections: https://github.com/broadinstitute/gnomad_qc/blob/main/gnomad_qc/example_notebooks/ancestry_classification_using_gnomad_rf.ipynb

import onnx
import hail as hl
from gnomad.sample_qc.ancestry import apply_onnx_classification_model, assign_population_pcs
from gnomad.utils.filtering import filter_to_adj

from gnomad_qc.v2.resources.basics import get_gnomad_meta
from gnomad_qc.v4.resources.basics import get_checkpoint_path

hl.init(driver_cores=8, worker_memory='highmem', tmp_dir="gs://schema_jsealock/tmp/")

MT = '/path/to/output/high_quality_sites.mt' # output from step 1
ANCESTRY_OUT = '/path/to/ancestry/output.ht'

# v3.1 PCA loadings and RF model.
gnomad_v3_loadings = "gs://gcp-public-data--gnomad/release/3.1/pca/gnomad.v3.1.pca_loadings.ht"
gnomad_v3_onnx_rf = "gs://gcp-public-data--gnomad/release/3.1/pca/gnomad.v3.1.RF_fit.onnx"

v3_num_pcs = 16
v3_min_prob = 0.75


with hl.hadoop_open(gnomad_v3_onnx_rf, "rb") as f:
    v3_onx_fit = onnx.load(f)


v3_loading_ht = hl.read_table(gnomad_v3_loadings)

# Filter to variants in the loadings Table.
mt = hl.read_matrix_table(MT)
mt = mt.semi_join_rows(v3_loading_ht)

# Project new genotypes onto loadings.
v3_pcs_ht = hl.experimental.pc_project(
    mt.GT, v3_loading_ht.loadings, v3_loading_ht.pca_af,
)

ht, model = assign_population_pcs(
    v3_pcs_ht,
    pc_cols=v3_pcs_ht.scores[:v3_num_pcs],
    fit=v3_onx_fit,
    min_prob=v3_min_prob,
    apply_model_func = apply_onnx_classification_model,
)

ht.write(ANCESTRY_OUT)
