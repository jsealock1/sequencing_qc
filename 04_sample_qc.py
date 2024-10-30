# run sample qc
# https://broadinstitute.github.io/gnomad_methods/api_reference/sample_qc/filtering.html#gnomad.sample_qc.filtering.compute_stratified_metrics_filter
import hail as hl
import gnomad
from gnomad.sample_qc.filtering import compute_stratified_metrics_filter
hl.init(driver_cores=8, worker_memory='highmem')

MT = '/path/to/output/variant/qced/matrix/table.mt' # output from step 1
ANCESTRY_OUT = '/path/to/ancestry/output.ht'
PASSING_SAMPLES = '/path/to/output/from/sample/qc.ht'

# read in variant qc'ed mt
mt = hl.read_matrix_table(MT)
ancestry = hl.read_tablee(ANCESTRY_OUT)

sample_qc_ht = hl.sample_qc(mt).cols().flatten().key_by('s')
ht = sample_qc_ht.annotate(qc_pop = ancestry[sample_qc_ht.s].Population)

filtering_qc_metrics = ['sample_qc.r_ti_tv','sample_qc.n_singleton','sample_qc.n_insertion','sample_qc.n_deletion','sample_qc.n_transition',
                        'sample_qc.n_transversion','sample_qc.r_het_hom_var','sample_qc.r_insertion_deletion']
stratified_metrics_ht = compute_stratified_metrics_filter(
                ht,
                qc_metrics={metric: ht[metric] for metric in filtering_qc_metrics},
                strata={"qc_pop": ht.qc_pop},
                metric_threshold={'sample_qc.n_singleton' : (4.0, 8.0)} ## we'll change the singleton upper bound to 8, to better fit a 0-bounded distribution
            )
passing_sample_qc = stratified_metrics_ht.filter((stratified_metrics_ht['fail_sample_qc.r_ti_tv']==False) & 
                                        (stratified_metrics_ht['fail_sample_qc.n_singleton']==False) & 
                                        (stratified_metrics_ht['fail_sample_qc.n_insertion']==False) & 
                                        (stratified_metrics_ht['fail_sample_qc.n_deletion']==False) & 
                                        (stratified_metrics_ht['fail_sample_qc.n_transition']==False) & 
                                        (stratified_metrics_ht['fail_sample_qc.n_transversion']==False) & 
                                        (stratified_metrics_ht['fail_sample_qc.r_het_hom_var']==False) & 
                                        (stratified_metrics_ht['fail_sample_qc.r_insertion_deletion']==False)
                                            )
mt = mt.filter_cols(hl.is_defined(passing_sample_qc[mt.col_key]))
passing_sample_qc.cols().write(PASSING_SAMPLES)
