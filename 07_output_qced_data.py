## output qc'ed matrix table
MT_PATH = '/path/to/matrix/table/of/called/genotypes.mt' # start with original matrix table
MT_OUT = '/path/to/output/with/variants/and/samples/qced.mt'

# passing samples from each step 
UNRELATED_SAMPLES = '/path/to/output/unrelated/samples.ht'
ANCESTRY_OUT = '/path/to/ancestry/output.ht'
PASSING_SAMPLES = '/path/to/output/from/sample/qc.ht'
SEX_CHECK_PASSING = '/path/to/output.ht'

# regions to filter (same as step 1)
TARGET_INTERVALS = "/target/intervals/path" 
LCR_PATH = "/lcr/path"

# VQSR results
VQSR_HT = '/path/to/variants/passing/vqsr'


# import original matrix table
mt = hl.read_matrix_table(MT_PATH)

# filter to passig samples
unrelated = hl.read_table(UNRELATED_SAMPLES)
ancestry = hl.read_table(ANCESTRY_OUT)
sample_qc = hl.read_table(PASSING_SAMPLES)
sex_check = hl.read_table(SEX_CHECK_PASSING)

mt = mt.filter_cols((hl.is_defined(unrelated[mt.col_key])) & (hl.is_defined(ancestry[mt.col_key])) & (hl.is_defined(sample_qc[mt.col_key])) & (hl.is_defined(sex_check[mt.col_key])))

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


# filter to sites passing VQSR
vqsr = hl.read_table(VQSR_HT) 
mt = mt.filter_rows(hl.is_defined(vqsr[mt.row_key]), keep=True)

mt_out.write(MT_OUT)
