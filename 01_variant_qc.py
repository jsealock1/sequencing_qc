## code for wes qc tutorial

import hail as hl
hl.init(driver_cores=8, worker_memory='highmem')

# inputs
MT_PATH = '/path/to/matrix/table/of/called/genotypes.mt'
TARGET_INTERVALS = "/target/intervals/path" 
LCR_PATH = "LCRFromHengHg38.bed"
VQSR_HT = '/path/to/variants/passing/vqsr'

# outputs
MT_OUT = '/path/to/output/variant/qced/matrix/table'

# import matrix table
mt = hl.read_matrix_table(MT_PATH)
# for this example we'll use the 1kgp + hgdp ref panel built into hail
mt = hl.experimental.load_dataset(name='gnomad_hgdp_1kg_subset_dense',
                                  version='3.1.2',
                                  reference_genome='GRCh38',
                                  region='us',
                                  cloud='gcp') 

# split multi allelic sites and add onto biallelic sites
bi = mt.filter_rows(hl.len(mt.alleles) == 2)
bi = bi.annotate_rows(a_index=1, was_split=False)
multi = mt.filter_rows(hl.len(mt.alleles) > 2)
split = hl.split_multi_hts(multi)
mt = split.union_rows(bi)

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

# make into minimal representation 
mt = mt.annotate_rows(min_rep = hl.min_rep(mt.locus, mt.alleles))
mt = mt.key_rows_by('min_rep')
mt = mt.drop('locus', 'alleles')

mt = mt.annotate_rows(locus = mt.min_rep.locus,
                      alleles = mt.min_rep.alleles)
mt = mt.key_rows_by('locus', 'alleles')
mt = mt.drop('min_rep')

# write MT
mt.write(MT_OUT)
