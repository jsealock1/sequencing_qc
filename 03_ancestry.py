# PCA

import hail as hl
hl.init(driver_cores=8, worker_memory='highmem')

# use gnomad projections: https://github.com/broadinstitute/gnomad_qc/blob/main/gnomad_qc/example_notebooks/ancestry_classification_using_gnomad_rf.ipynb

MT = '/path/to/output/variant/qced/matrix/table.mt' # output from step 1
META = "/path/to/case/control/status"
PRUNED_VARIANT_TABLE = "/path/to/ld/pruned/variants.ht"
PASSING_SAMPLES = '/path/to/unrelated/samples.ht'
ANCESTRY_OUT = '/path/to/ancestry/output.ht'
PCA_OUT = '/path/to/write/pc/output.ht'

# read in variant qc'ed mt
mt = hl.read_matrix_table(MT)

## filter to unrelated samples
rel_samples_keep = hl.read_table(PASSING_SAMPLES)
mt = mt.filter_cols(hl.is_defined(rel_samples_keep[mt.col_key]))

##load ref data
gnomad1 = hl.experimental.load_dataset(name='gnomad_hgdp_1kg_subset_dense',
                                  version='3.1.2',
                                  reference_genome='GRCh38',
                                  region='us',
                                  cloud='gcp') # https://hail.is/docs/0.2/experimental/index.html#hail.experimental.load_dataset
gnomad = gnomad1.annotate_cols(Dataset = 'gnomad_hgdp_1kg')
# Take subset of entry fields 
gnomad = gnomad.select_entries('DP', 'GQ', 'GT')
# Take only sites passing VQSR
gnomad = gnomad.filter_rows(gnomad.filters == hl.empty_set(hl.tstr), keep = True)
# Drop other column/global annotations
gnomad = gnomad.select_cols(gnomad.Dataset)
gnomad = gnomad.select_rows('rsid')
gnomad = gnomad.select_globals()

# combine datasets
pca_mt = mt.union_cols(gnomad)

# filter to ld pruned variants
pruned_variant_table = hl.read_table(PRUNED_VARIANT_TABLE)
pca_mt = mt.filter_rows(hl.is_defined(pruned_variant_table[pca_mt.row_key]))

# remove variants with high missingness
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.call_rate>0.99)

# run pca
eigenvalues, score_table, loading_table = hl.hwe_normalized_pca(pca_mt.GT, k=10, compute_loadings=True)
score_table.write(PCA_OUT)

# random forest for ancestry 
# https://broadinstitute.github.io/gnomad_methods/api_reference/sample_qc/ancestry.html#gnomad.sample_qc.ancestry.assign_population_pcs
