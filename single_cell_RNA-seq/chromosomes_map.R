# chromoMap provides interactive, configurable and elegant graphics visualization of chromosomes or chromosomal regions allowing users to map chromosome features (like genes,SNPs etc.) onto the chromosome and visualize the feature-associated data (like multi-omics data). Each chromosome is composed of genomic-windows(representing a specific range determined based on chromosome length) that, on hover, shows details about the annotations in that locus range. The plots can be saved as HTML documents that can be shared easily. In addition, you can include them in R Markdown or in R Shiny applications
# https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html#Session_Info
##### Install chromoMap #####
install.packages("chromoMap")
install.packages("RIdeogram")
??chromosome_map
##### Prepare Input Data Files #####
# chromoMap files for this vignette
# chromosome files
rm(list = ls())
getwd()
chr_file_1 = "./single_cell_RNA-seq/chr_file_without_centromere.txt"
chr_file_2 = "./single_cell_RNA-seq/chr_file_with_centromere.txt"


# annotation files
anno_file_1 <- "./single_cell_RNA-seq/annotation_pos.txt"
anno_file_2 <- "./single_cell_RNA-seq/annotation_pos_and_neg.txt"

df = data.frame(chr=c("chr16"),
                start=2475127,end=c(2505730),cent=c(2495030))
write.table(df,"C:\\Users\\wane\\Documents\\file\\sci\\Shuntak\\MR\\chromosome_file.txt",quote=F,
            sep="\t",row.names=FALSE,col.names=FALSE)
anno = data.frame(Ele=c("An1","An2"),
                  chr=c("chr16","chr16"),start=c(2485127,2551127),
                  end=c(2485573,2551730),score=c(10,177))
write.table(anno,"C:\\Users\\wane\\Documents\\file\\sci\\Shuntak\\MR\\anno.txt",quote=F,
            sep="\t",row.names=FALSE,col.names=FALSE)

library(chromoMap)
chromoMap("C:\\Users\\wane\\Documents\\file\\sci\\Shuntak\\MR\\chromosome_file.txt",
          "C:\\Users\\wane\\Documents\\file\\sci\\Shuntak\\MR\\anno.txt")

# R objects as chromoMap Inputs
# passing data.frames directly instead of files
chromoMap(list(chr.data), list(anno.data))
# for polyploidy
chromoMap(list(chr.data1, chr.data2),
  list(anno.data1, anno.data2),
  ploidy = 2
)

# My first chromoMap
chromoMap(chr_file_1, anno_file_1)

# Now, let’s create one with centromeres.
chromoMap(chr_file_2,anno_file_1)

# Point and Segment-annotation plots
chromoMap("./single_cell_RNA-seq/chr_file_with_centromere.txt",
          "./single_cell_RNA-seq/annotation_pos_and_neg.txt",
          segment_annotation = T)






























