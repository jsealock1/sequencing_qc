import hail as hl
hl.init(driver_cores=8, worker_memory='highmem')

MT = '/path/to/output/variant/qced/matrix/table.mt' # output from step 1
META = '/path/to/meta/data/with/gender'
SEX_CHECK_PASSING = '/path/to/output.ht'

# read in variant qc'ed mt
mt = hl.read_matrix_table(MT)

# filter to bialleleic chr x variants
mt_x = mt.filter_rows(mt.locus.contig == 'chrX') 
dataset = mt_x.filter_rows(hl.len(mt.alleles) == 2)
dataset = hl.variant_qc(dataset)
dataset = dataset.filter_rows(dataset.qc.AF[1] >= 0.01)

# ld prune
pruned_variant_table = hl.ld_prune(dataset.GT)
dataset = dataset.filter_rows(hl.is_defined(pruned_variant_table[dataset.row_key]))

# impute sex
imputed_sex = hl.impute_sex(dataset.GT, aaf_threshold=0.01, include_par=False) 

# check with reported gender
imputed_sex = imputed_sex.annotate(reported_gender = meta[imputed_sex.s].gender)
pass_sex_check = imputed_sex.filter(((imputed_sex.reported_gender=="M") & (imputed_sex.impute_sex.is_female=="false")) | ((imputed_sex.reported_gender=="F") & (imputed_sex.impute_sex.is_female=="true")))

pass_sex_check.write(SEX_CHECK_PASSING)