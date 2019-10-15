#!/usr/bin/Rscript

# imports ==================
library(VariantAnnotation)
library(readr)

# params =====================
args = commandArgs(trailingOnly=TRUE)
vcf_file = args[1]
cna_fn = args[2]
tmp_file = args[3]
cna_type = args[4]
# this is the sample in the vcf. some callers have TUMOUR, others other
# names
vcf_sample = args[5]

prefix = stringr::str_extract(string = tmp_file,
                              pattern = "[:graph:]+(?=.tmpFile)")

print(prefix)




# function parsers --------------------------------------------------------

read_hartwig <- function(cna_fn){
  dat = readr::read_tsv(cna_fn, col_types = cols(
    `#chromosome` = col_character(),
    start = col_double(),
    end = col_double(),
    copyNumber = col_double(),
    bafCount = col_double(),
    observedBAF = col_double(),
    baf = col_double(),
    segmentStartSupport = col_character(),
    segmentEndSupport = col_character(),
    method = col_character(),
    depthWindowCount = col_double(),
    gcContent = col_double(),
    minStart = col_double(),
    maxStart = col_double(),
    minorAllelePloidy = col_double(),
    majorAllelePloidy = col_double()
  ))

  res_dat = data.frame(
    chromosome = dat$`#chromosome`,
    start = dat$start,
    end = dat$end,
    major_cn = round(dat$majorAllelePloidy),
    minor_cn = round(dat$minorAllelePloidy)

  )
  return(res_dat)

}

read_sanger_pcawg <- function(cna_fn) {
  cna_dat = readr::read_tsv(cna_fn, col_types =cols(
    chromosome = col_character(),
    start = col_double(),
    end = col_double(),
    total_cn = col_double(),
    major_cn = col_double(),
    minor_cn = col_double(),
    star = col_double()
  )
  )

  return(cna_dat)
}

# script ------------------------------------------------------------------

cna_parsers = c("sanger_pcawg","hartwig")

stopifnot(cna_type %in% cna_parsers)

cna_dat = switch (cna_type,
                  sanger_pcawg  = read_sanger_pcawg(cna_fn = cna_fn),
                  hartwig = read_hartwig(cna_fn = cna_fn)
)

# some CNA have NA values
mask_na = is.na(cna_dat$major_cn) | is.na(cna_dat$major_cn)
cna_dat = cna_dat[!mask_na,]

# pyClone also doesn't like CPmajor = 0, this should be a pretty rare
# case so I also remove it.
mask_cp0 =  cna_dat$major_cn == 0
cna_dat = cna_dat[!mask_cp0,]

if (any(mask_cp0)){
  warning(glue::glue("Some regions removed due to CNA Major = 0"))
}

cna_gr = GRanges(seqnames = cna_dat$chromosome,
                  ranges = IRanges(start = cna_dat$start,
                   end = cna_dat$end) )

vcf_dat = VariantAnnotation::readVcf(vcf_file)
vcf_vr = as(object= vcf_dat, Class= "VRanges" )

# remove normal sample
vcf_vr = vcf_vr[sampleNames(vcf_vr) == vcf_sample]

# remove sex chromosomes (?) =============
vcf_vr = vcf_vr[seqnames(vcf_vr) %in%  1:22]

# overlap with the CNA
ovr = findOverlaps(query = vcf_vr,subject = cna_gr)
# remove mutations with no available overlap
# I check and they seem to be mutations that outside of the
# CNA called
vcf_vr = vcf_vr[queryHits(ovr)]

# transform to pyClone format =============
pyclone_df = data.frame(
    mutation_id = glue::glue(
      "{seqnames(vcf_vr)}:{end(vcf_vr)}:{ref(vcf_vr)}:{alt(vcf_vr)}"),
    minor_cn = cna_dat[["minor_cn"]][subjectHits(ovr)],
    major_cn = cna_dat[["major_cn"]][subjectHits(ovr)])


preproc_file = readr::read_tsv(tmp_file)

res_df = dplyr::left_join(pyclone_df,preproc_file)

# reorder the columns of the resulting data frame
res_df = res_df[,c("mutation_id",
                   "ref_counts",
                   "var_counts",
                   "normal_cn",
                   "minor_cn",
                   "major_cn")]

# mutation_id
# ref_counts
# var_counts
# normal_cn
# minor_cn
# major_cn

readr::write_tsv(x= res_df,
                 path= glue::glue("{prefix}_pyclone_format.tsv"))
