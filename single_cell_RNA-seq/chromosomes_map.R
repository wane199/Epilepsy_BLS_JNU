# chromoMap provides interactive, configurable and elegant graphics visualization of chromosomes or chromosomal regions allowing users to map chromosome features (like genes,SNPs etc.) onto the chromosome and visualize the feature-associated data (like multi-omics data). Each chromosome is composed of genomic-windows(representing a specific range determined based on chromosome length) that, on hover, shows details about the annotations in that locus range. The plots can be saved as HTML documents that can be shared easily. In addition, you can include them in R Markdown or in R Shiny applications
# https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html#Session_Info
##### Install chromoMap #####
install.packages("chromoMap")
install.packages("RIdeogram")

##### Prepare Input Data Files #####
# chromoMap files for this vignette
# chromosome files
# chr_file_1 = "chr_file_without_centromere.txt"
# chr_file_2 = "chr_file_with_centromere.txt"
getwd()
chr_file_1 <- read.table("C:\\Users\\wane\\Documents\\file\\sci\\Shuntak\\MR\\chr_file_without_centromere.csv", sep = ";")
chr_file_2 <- read.table("C:\\Users\\wane\\Documents\\file\\sci\\Shuntak\\MR\\chr_file_with_centromere.csv", sep = ";")

# annotation files
anno_file_1 <- read.table("C:\\Users\\wane\\Documents\\file\\sci\\Shuntak\\MR\\annotation_pos.csv", sep = ";")
anno_file_2 <- read.table("C:\\Users\\wane\\Documents\\file\\sci\\Shuntak\\MR\\annotation_pos_and_neg.csv", sep = ";")

# R objects as chromoMap Inputs
# passing data.frames directly instead of files
chromoMap(list(chr.data), list(anno.data))
# for polyploidy
chromoMap(list(chr.data1, chr.data2),
  list(anno.data1, anno.data2),
  ploidy = 2
)

# My first chromoMap
library(chromoMap)
chromoMap(chr_file_1, anno_file_1)

# Now, let’s create one with centromeres.
chromoMap(chr_file_2,anno_file_1)































