#!/usr/bin/Rscript

# David Mas
# Compute clonality of mutations from pyclone format

# functions ---------------------------------------------------------------
compute_clonality_pyClone_format <- function(muts_df,
                                             purity,
                                             CPnorm){
  ## check colnames
  req_cols = c("mutation_id",
               "ref_counts",
               "var_counts",
               "normal_cn",
               "minor_cn",
               "major_cn")

  stopifnot(all(colnames(muts_df) == req_cols))

  ## compute expected VAF
  qt = (muts_df$major_cn + muts_df$minor_cn)
  CPnorm = 2

  Sq = apply(matrix(c(muts_df$major_cn,
                      muts_df$minor_cn),
                    byrow = FALSE,
                    ncol = 2),
             1,
             max)
  #browser()
  Sq[Sq == 0] = 1

  muts_df$EVaf = (Sq*purity) / (2 * (1 - purity) + purity*qt)

  # for some reason this happens...
  impossible_mask = (muts_df$major_cn == muts_df$minor_cn & muts_df$minor_cn == 0)
  if (sum(impossible_mask) > 0){
    warning(glue::glue("IMP events {sum(impossible_mask)}, setting to 0.5"))
    muts_df[impossible_mask]$EVaf = 0.5 # (?)
  }

  # perform the test.
  res = binom_test(x = muts_df$var_counts,
                   n = (muts_df$ref_counts + muts_df$var_counts),
                   p = muts_df$EVaf,
                   alternative = "less")
  res$qval = p.adjust(p = res$p.value,method = "fdr")
  muts_df = cbind(muts_df,res)

  return(muts_df)

}

# params ========
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-i", "--inputFile"),
    action = "store",
    default = NA,
    type = 'character',
    help = "a file in pyClone format"
  ),
  make_option(
    c("-p", "--purity"),
    action = "store",
    default = NA,
    type = 'double',
    help = "the cancer cell fraction of the sample"
  ),
  make_option(
    c("-v", "--verbose"),
    action = "store_true",
    default = TRUE,
    help = "Should the program print extra stuff out? [default %default]"
  ),
  make_option(
    c("-q", "--quiet"),
    action = "store_false",
    dest = "verbose",
    help = "Make the program not be verbose."
  ),
  make_option(
    c("-P", "--plot"),
    action = "store_true",
    default = FALSE,
    help = "Flag to plot results [default %default]"
  ),
  make_option(
    c("-s", "--sampleName"),
    action = "store",
    default = "output",
    type = 'character',
    help = "Prefix for generated outputs. [default %default]"
  )
)
opt = parse_args(OptionParser(option_list=option_list))


## todo
##
## add verbose output
## clean output in directory

## here for dbug
if (interactive()){
  opt$inputFile = "Y:/users/dmas/data/PCAWG_MUTS/SANGER_WGS/pyClone/tmp_WD/42/92a6b8353098ed17732d02de445657/fc9404ed-1ba3-2638-e040-11ac0c484da2_pyclone_format.tsv"
  opt$purity = 0.8
  opt$plot = TRUE
  opt$sampleName = "fc9404ed-1ba3-2638-e040-11ac0c484da2"
}


# imports ---------------
library(magrittr)
library(ggplot2)
library(VariantAnnotation)
library(genomicHelpersDMP)



# data --------------------------------------------------------------------

# check if files are available and exists
if (!fs::file_exists(opt$inputFile)){
  stop("No files provided.")
}

muts_df = readr::read_tsv(file = opt$inputFile)

res_df = compute_clonality_pyClone_format(muts_df = muts_df,
                                          purity = opt$purity,
                                          CPnorm = 2)

# estimate groups with mclust


# script ------------------------------------------------------------------

library(mclust)

ratio = log(res_df$estimate / res_df$EVaf)

ol = length(ratio)
incompatible_mask = is.na(ratio) | is.infinite(ratio)
ratio = ratio[!incompatible_mask]
res_df = res_df[!incompatible_mask,]

per = scales::percent(sum(incompatible_mask)/ol)
warning(glue::glue("{per} mutations were incompatible, removing those."))

fs::dir_create("clonality_plots")
set.seed(42)
res_df$ratio = ratio

#add opt
sname = opt$sampleName

model_gs = mclustBIC(data = res_df$ratio,
                         G = 1:4,
                         modelNames = c("E"))
clonality = Mclust(data = res_df$ratio,
                       x = model_gs)

res_df$classification = clonality$classification



# OUTPUT ------------------------------------------------------------------
print("output start")

plot_df = data.frame(ratio = res_df$ratio,type = factor(res_df$classification))
print(nrow(plot_df))
plot = plot_df %>%
    ggplot(aes(x = ratio,fill = factor(type))) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = 0,linetype = "dashed") +
    labs(title = sname,
         x = "ln(obs/exp)") +
  theme_classic() +
  scale_color_brewer(palette = "Set1")

ggsave(filename = glue::glue("clonality_plots/{sname}_clonality_ratio.pdf"),
         plot = plot,
         device = "pdf")

# save results ==================
fs::dir_create("clonality_results")

# we need this to repaste the sample name in the middle of the next 
# id that is recognised by clustmut. This is a small design flaw. 
# maybe I should open an issue
res_df$mutation_id %>% stringr::str_split(":") -> mId_list
seqname_vec = mId_list %>% purrr::map_chr(1)
start_vec = mId_list %>% purrr::map_chr(2)
ref_vec = mId_list %>% purrr::map_chr(3)
alt_vec = mId_list %>% purrr::map_chr(4)

mid_vec = glue::glue("{seqname_vec}:{start_vec}:{sname}:{ref_vec}:{alt_vec}")
y = res_df

clonality = glue::glue("clonality{y$classification}")

# pair + clonality
pair = ifelse(ref_vec %in% c("C","G"),yes = "S",no = "W")
string_code = glue::glue("{mid_vec}_{clonality}{pair}")
readr::write_lines(string_code,
                       path = glue::glue("clonality_results/{sname}_mutations_pair_clonality.txt"))

# pair + strand
string_code = glue::glue("{mid_vec}_{clonality}{ref_vec}")
readr::write_lines(string_code,
                       path = glue::glue("clonality_results/{sname}_mutations_strand_clonality.txt"))

# only clonality
string_code = glue::glue("{mid_vec}_{clonality}")
    readr::write_lines(string_code,
                       path = glue::glue("clonality_results/{sname}_mutations_clonality.txt"))

