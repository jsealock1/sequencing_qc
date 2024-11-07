import hail as hl
hl.init(driver_cores=8, worker_memory='highmem')

MT = '/path/to/output/high_quality_sites.mt' # output from step 1
META = '/path/to/meta/data/with/gender'
SEX_CHECK_PASSING = '/path/to/output.ht'

# read in variant qc'ed mt
mt = hl.read_matrix_table(MT)

# filter to bialleleic chr x variants
mt_x = mt.filter_rows(mt.locus.contig == 'chrX') 
mt_x = mt_x.filter_rows(hl.len(mt.alleles) == 2)
mt_x = hl.variant_qc(mt_x)
mt_x = mt_x.filter_rows(mt_x.qc.AF[1] >= 0.01)

# impute sex
imputed_sex = hl.impute_sex(mt_x.GT, aaf_threshold=0.01, include_par=False) 

# check with reported gender
imputed_sex = imputed_sex.annotate(reported_gender = meta[imputed_sex.s].gender)
pass_sex_check = imputed_sex.filter(((imputed_sex.reported_gender=="M") & (imputed_sex.impute_sex.is_female=="false")) | ((imputed_sex.reported_gender=="F") & (imputed_sex.impute_sex.is_female=="true")))

pass_sex_check.write(SEX_CHECK_PASSING)