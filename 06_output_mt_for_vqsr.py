## output qc'ed matrix table
MT_PATH = '/path/to/matrix/table/of/called/genotypes.mt' # start with original matrix table

# passing samples from each step 
UNRELATED_SAMPLES = '/path/to/output/unrelated/samples.ht'
ANCESTRY_OUT = '/path/to/ancestry/output.ht'
PASSING_SAMPLES = '/path/to/output/from/sample/qc.ht'
SEX_CHECK_PASSING = '/path/to/output.ht'

# output
MT_FOR_VQSR = '/path/to/matrix/table/for/vqsr.mt'


# import matrix table
mt = hl.read_matrix_table(MT_PATH)

unrelated = hl.read_table(UNRELATED_SAMPLES)
ancestry = hl.read_table(ANCESTRY_OUT)
sample_qc = hl.read_table(PASSING_SAMPLES)
sex_check = hl.read_table(SEX_CHECK_PASSING)

mt = mt.filter_cols((hl.is_defined(unrelated[mt.col_key])) & (hl.is_defined(ancestry[mt.col_key])) & (hl.is_defined(sample_qc[mt.col_key])) & (hl.is_defined(sex_check[mt.col_key])))
mt.write(MT_FOR_VQSR)
## export this matrix table and run vqsr (https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR)

