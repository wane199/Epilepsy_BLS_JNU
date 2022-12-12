## Total volumes of distribution (VT) from a simulated PET study
## including a baseline scan, as well as two other scans after
## administration of a drug. Note that each row in the matrix
## represents a ROI, whilst each column represents a scan.
library(occ)
data(occ.example)
occ.example

m = occ(occ.example)

print(m)     # Prints the neuroreceptor occupancy coefficients

summary(m)   # Also prints the non-displaceable volume of

fitted(m)    # Prints the fitted values

residuals(m) # Prints the residuals

plot(m)      # Plots the estimated and observed volumes of
# distribution


library(oro.pet)
n <- 11
h <- seq(200, 150, length=n)
w <- seq(80, 120, length=n)
cbind(h, w, leanBodyMass(h, w, "male"), leanBodyMass(h, w, "female"))


leanBodyMass(height, weight, gender)

hotSpotSUV(suv, radius = 10, type = "3D")

totalSUV(suv, mask, z, bg, local = TRUE)


#############################
# wordcloud
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) # 设置清华镜像
install.packages("jsonlite")    #一定要安装这个包
library('devtools')
devtools::install_github("lchiffon/wordcloud2")

# library
library(wordcloud2) 
letterCloud( demoFreq, word = "THANKS PET!", color='random-light' , backgroundColor="black")
my_graph <- letterCloud( demoFreq, word = "PET", color="white", backgroundColor="pink")

# Export the wordcloud
# Wordcloud2 is a html widget. It means your wordcloud will be output in a HTML format.You can export it as a png image using rstudio, or using the webshot library as follow:

# load wordcloud2
library(wordcloud2) 
# install webshot
library(webshot)
webshot::install_phantomjs()

# Make the graph
my_graph <- wordcloud2(demoFreq, size=1.5)

# save it in html
library("htmlwidgets")
saveWidget(my_graph,"./BLS_EP_files/tmp.html",selfcontained = F)

# and in png or pdf
webshot("./BLS_EP_files/tmp.html","./BLS_EP_files/fig_1.pdf", delay =5, vwidth = 480, vheight=480)
