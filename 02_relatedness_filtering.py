
import hail as hl
hl.init(driver_cores=8, worker_memory='highmem')

MT = '/path/to/output/high_quality_sites.mt' # output from step 1
META = "/path/to/case/control/status"
PRUNED_VARIANT_TABLE = "/path/to/output/ld/pruned/variants.ht"
IBD_OUTPUT = '/path/to/output/from/relatedness.ht'
UNRELATED_SAMPLES = '/path/to/output/unrelated/samples.ht'

# read in variant qc'ed mt
mt = hl.read_matrix_table(MT)

# relatedness filtering 
meta = hl.import_table(META, key='s')
mt = mt.annotate_cols(status = meta[mt.s].status)
mt = mt.filter_cols((mt.status=="case") | (mt.status=="control"), keep=True)

# run pc_relate
mt = mt.annotate_cols(is_case = hl.if_else(mt.status=="case", hl.bool(True), hl.bool(False)))
pc_rel = hl.pc_relate(mt.GT, 0.01, k=10, statistics='kin', min_kinship=0.2)
pairs = pc_rel.filter(pc_rel['kin'] > 0.2)

# filter for max independent set for case/control
# https://hail.is/docs/0.2/methods/misc.html#hail.methods.maximal_independent_set
samples = mt.cols()
pairs_with_case = pairs.key_by(
    i=hl.struct(id=pairs.i, is_case=samples[pairs.i].is_case),
    j=hl.struct(id=pairs.j, is_case=samples[pairs.j].is_case))
def tie_breaker(l, r):
    return hl.if_else(l.is_case & ~r.is_case, -1,
                      hl.if_else(~l.is_case & r.is_case, 1, 0))

related_samples_to_remove = hl.maximal_independent_set(pairs_with_case.i, pairs_with_case.j, False, tie_breaker)
result = mt.filter_cols(hl.is_defined(related_samples_to_remove.key_by(s = related_samples_to_remove.node.id.s)[mt.col_key]), keep=False)
rel_samples_keep = result.cols()

# save data
pairs.write(IBD_OUTPUT)
rel_samples_keep.write(UNRELATED_SAMPLES)
