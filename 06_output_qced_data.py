## output qc'ed matrix table
MT = '/path/to/output/variant/qced/matrix/table.mt' # output from step 1
MT_OUT = '/path/to/output/with/variants/and/samples/qced.mt'
UNRELATED_SAMPLES = '/path/to/output/unrelated/samples.ht'
ANCESTRY_OUT = '/path/to/ancestry/output.ht'
PASSING_SAMPLES = '/path/to/output/from/sample/qc.ht'
SEX_CHECK_PASSING = '/path/to/output.ht'

mt = hl.read_matrix_table(MT)

unrelated = hl.read_table(UNRELATED_SAMPLES)
ancestry = hl.read_table(ANCESTRY_OUT)
sample_qc = hl.read_table(PASSING_SAMPLES)
sex_check = hl.read_table(SEX_CHECK_PASSING)

mt_out = mt.filter_cols((hl.is_defined(unrelated[mt.col_key])) & (hl.is_defined(ancestry[mt.col_key])) & (hl.is_defined(sample_qc[mt.col_key])) & (hl.is_defined(sex_check[mt.col_key])))
mt_out.write(MT_OUT)
