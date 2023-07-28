## code for wes qc tutorial

import hail as hl

hl.init(driver_cores=8, worker_memory='highmem', tmp_dir="gs://schema_jsealock/tmp/")

# inputs
MT_PATH = '/path/to/matrix/table/of/called/genotypes.mt'
TARGET_INTERVALS = "Twist_Alliance_Clinical_Research_Exome_Covered_Targets_hg38-34.9MB.bed" 
LCR_PATH = "LCRFromHengHg38.bed"
VQSR_HT = '/path/to/variants/passing/vqsr'
PCA_SNPS = "purcell5k_grch38_liftover_2021-09-14.interval_list"
META = "/path/to/case/control/status"

# outputs
MT_OUT = '/path/to/output/qced/matrix/table'
IBD_OUTPUT = '/path/to/output/from/relatedness'
PASSING_SAMPLES = '/path/to/output/from/sample/qc'

# import matrix table
mt = hl.read_matrix_table(MT_PATH)
# for this example we'll use the 1kgp + hgdp ref panel built into hail
mt = hl.experimental.load_dataset(name='gnomad_hgdp_1kg_subset_dense',
                                  version='3.1.2',
                                  reference_genome='GRCh38',
                                  region='us',
                                  cloud='gcp') 


# VQSR filter
## in this dataset, VQSR results are stored within the matrix table 
mt = mt.filter_rows(mt.filters == hl.empty_set(hl.tstr), keep = True)

## however, most experimental datasets will have a separate file with VQSR results to be loaded and filtered:
vqsr = hl.read_table(VQSR_HT) 
mt = mt.filter_rows(hl.is_defined(vqsr[mt.row_key]), keep=True)

#filter out LCRs 
lcr_intervals = hl.import_locus_intervals(LCR_PATH, reference_genome='GRCh38', skip_invalid_intervals=True)
mt = mt.filter_rows(hl.is_defined(lcr_intervals[mt.locus]), keep=False)

# interval filter for exome
intervals = hl.import_locus_intervals(TARGET_INTERVALS, reference_genome="GRCh38")
mt = mt.annotate_rows(not_in_target_intervals = ~hl.is_defined(intervals[mt.locus]))
mt = mt.filter_rows(mt.not_in_target_intervals, keep=False)


# variant QC
# filter for GQ, DP, and AB
mt = mt.filter_entries(
    hl.is_defined(mt.GT) &
    (
        (mt.GT.is_hom_ref() & 
            (
                (mt.GQ < 20) |
                (mt.DP < 10)
        	)
        ) |
        (mt.GT.is_het() & 
        	( 
                (((mt.AD[0] + mt.AD[1]) / mt.DP) < 0.8) | 
                ((mt.AD[1] / mt.DP) < 0.2) | 
                (mt.PL[0] < 20) |
                (mt.DP < 10)
        	)
        ) |
        (mt.GT.is_hom_var() & 
        	(
                ((mt.AD[1] / mt.DP) < 0.8) |
                (mt.PL[0] < 20) |
                (mt.DP < 10)
        	)
        )
    ),
    keep = False
)

# Filter out the invariant rows.
mt = hl.variant_qc(mt, name='qc')
mt = mt.filter_rows((mt.qc.AF[0] > 0.0) & (mt.qc.AF[0] < 1.0))

# make into minimal representation 
mt = mt.annotate_rows(min_rep = hl.min_rep(mt.locus, mt.alleles))
mt = mt.key_rows_by('min_rep')
mt = mt.drop('locus', 'alleles')

mt = mt.annotate_rows(locus = mt.min_rep.locus,
                      alleles = mt.min_rep.alleles)
mt = mt.key_rows_by('locus', 'alleles')
mt = mt.drop('min_rep')

# relatedness
meta = hl.import_table(META, key='s')
mt = mt.annotate_cols(status = meta[mt.s].status)
mt = mt.filter_cols((mt.status=="case") | (mt.status=="control"), keep=True)

pca_snps = hl.import_locus_intervals(PCA_SNPS, reference_genome="GRCh38")
mt = mt.filter_rows(hl.is_defined(pca_snps[mt.locus]), keep=True)
mt = mt.select_rows('rsid')
mt = mt.select_entries('DP', 'GQ', 'GT')

dataset = mt
dataset = dataset.annotate_cols(is_case = hl.if_else(dataset.status=="case", hl.bool(True), hl.bool(False)))
pc_rel = hl.pc_relate(dataset.GT, 0.01, k=10, statistics='kin', min_kinship=0.2, block_size=2048)
pairs = pc_rel.filter(pc_rel['kin'] > 0.2)
pairs.write(IBD_OUTPUT)

# filter for max independent set for case/control
samples = dataset.cols()
pairs_with_case = pairs.key_by(
    i=hl.struct(id=pairs.i, is_case=samples[pairs.i].is_case),
    j=hl.struct(id=pairs.j, is_case=samples[pairs.j].is_case))
def tie_breaker(l, r):
    return hl.if_else(l.is_case & ~r.is_case, -1,
                      hl.if_else(~l.is_case & r.is_case, 1, 0))


related_samples_to_remove = hl.maximal_independent_set(
   pairs_with_case.i, pairs_with_case.j, False, tie_breaker)
result = dataset.filter_cols(hl.is_defined(
	related_samples_to_remove.key_by(
		s = related_samples_to_remove.node.id.s)[dataset.col_key]), keep=False)

samples_keep = result.cols()
samples_keep.export(PASSING_SAMPLES)

# PCA
## pca with passing samples
passing = hl.import_table(PASSING_SAMPLES, key="s")
mt = mt.filter_cols(hl.is_defined(passing[mt.col_key]))

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

# take purcell 5k snps from both datasets 
gnomad = gnomad.filter_rows(hl.is_defined(pca_snps[gnomad.locus]), keep=True)

# combine datasets
mt = mt.union_cols(gnomad)

pop = gnomad1.select_cols(gnomad1.gnomad_population_inference.pop, gnomad1.hgdp_tgp_meta.project).cols()
pop = pop.select_globals()
pop = pop.filter(pop.project=='synthetic_diploid_truth_sample', keep=False)
pop.flatten().export(f'{PCA_OUT}_gnomad_hgdp_1kgp_pop_labels.txt')

# run pca
eigenvalues, score_table, loading_table = hl.hwe_normalized_pca(mt.GT, k=20, compute_loadings=True)
score_table.flatten().export(f'{PCA_OUT}_pc_scores.txt')

# random forest for ancestry 

# sex imputation
### per ancestry impute_sex
### keep samples with same sex and reported gender

# run sample qc
mt = mt.filter_cols(hl.is_defined(ibd_passing[mt.col_key]))
mt = mt.filter_cols(hl.is_defined(pca_passing[mt.col_key]))

sample_qc = hl.sample_qc(mt)
sample_qc.cols().select('sample_qc').flatten().export(output=SAMPLE_QC_OUT)


# sample qc filtering
### per ancestry
### remove samples +/- MADs from median 


