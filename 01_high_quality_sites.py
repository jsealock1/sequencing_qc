import hail as hl
hl.init(driver_cores=8, worker_memory='highmem')

# inputs
MT_PATH = '/path/to/matrix/table/of/called/genotypes.mt'
TARGET_INTERVALS = "/target/intervals/path" 
LCR_PATH = "/lcr/path"

# outputs
MT_OUT = '/path/to/output/high_quality_sites.mt'

# import matrix table
mt = hl.read_matrix_table(MT_PATH)

#filter out LCRs 
lcr_intervals = hl.import_locus_intervals(LCR_PATH, reference_genome='GRCh38', skip_invalid_intervals=True)
mt = mt.filter_rows(hl.is_defined(lcr_intervals[mt.locus]), keep=False)

# interval filter for exome
intervals = hl.import_locus_intervals(TARGET_INTERVALS, reference_genome="GRCh38")
mt = mt.annotate_rows(not_in_target_intervals = ~hl.is_defined(intervals[mt.locus]))
mt = mt.filter_rows(mt.not_in_target_intervals, keep=False)

# filter for GQ, DP, and AB
mt = mt.filter_entries(
    hl.is_defined(mt.GT) &
    ((mt.GT.is_hom_ref() & ((mt.GQ < 20) | (mt.DP < 10))) |
    (mt.GT.is_het() & ((((mt.AD[0] + mt.AD[1]) / mt.DP) < 0.8) | ((mt.AD[1] / mt.DP) < 0.2) | (mt.PL[0] < 20) |(mt.DP < 10))) |
    (mt.GT.is_hom_var() & (((mt.AD[1] / mt.DP) < 0.8) | (mt.PL[0] < 20) | (mt.DP < 10)))),
    keep = False)

# Filter out the invariant rows.
mt = hl.variant_qc(mt, name='qc')
mt = mt.filter_rows((mt.qc.AF[0] > 0.0) & (mt.qc.AF[0] < 1.0))

## find common, highly called sites
mt = hl.variant_qc(mt)
mt = mt.filter_rows((mt.variant_qc.AF[1] >= 0.01) & (mt.variant_qc.call_rate >= 0.95))

# LD prune to independent sites
pruned_variants = hl.ld_prune(mt)
mt = mt.filter_rows(hl.is_defined(pruned_variants[mt.row_key]), keep=True)

# write MT
mt.write(MT_OUT)
