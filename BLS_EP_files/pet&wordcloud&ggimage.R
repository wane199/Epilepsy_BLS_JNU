## Total volumes of distribution (VT) from a simulated PET study
## including a baseline scan, as well as two other scans after
## administration of a drug. Note that each row in the matrix
## represents a ROI, whilst each column represents a scan.
library(occ)
data(occ.example)
occ.example

m <- occ(occ.example)

print(m) # Prints the neuroreceptor occupancy coefficients

summary(m) # Also prints the non-displaceable volume of

fitted(m) # Prints the fitted values

residuals(m) # Prints the residuals

plot(m) # Plots the estimated and observed volumes of
# distribution


library(oro.pet)
n <- 11
h <- seq(200, 150, length = n)
w <- seq(80, 120, length = n)
cbind(h, w, leanBodyMass(h, w, "male"), leanBodyMass(h, w, "female"))


leanBodyMass(height, weight, gender)

hotSpotSUV(suv, radius = 10, type = "3D")

totalSUV(suv, mask, z, bg, local = TRUE)


#############################
# wordcloud
options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) # 设置清华镜像
install.packages("jsonlite") # 一定要安装这个包
library("devtools")
devtools::install_github("lchiffon/wordcloud2")

# library
library(wordcloud2)
head(demoFreq)
letterCloud(ai,
  word = "PET",
  color = "random-light",
  backgroundColor = "black"
)
ai <- read.csv("./BLS_EP_files/Most_Frequent_Words.csv")
rownames(ai) <- ai[, 1]
head(ai)
letterCloud(ai,
  word = "PET",
  color = "random-light",
  backgroundColor = "pink"
)
my_graph <-
  letterCloud(
    ai,
    word = "PET",
    color = "white",
    backgroundColor = "pink"
  )
my_graph

# Export the wordcloud
# Wordcloud2 is a html widget. It means your wordcloud will be output in a HTML format.You can export it as a png image using rstudio, or using the webshot library as follow:
# load wordcloud2
library(wordcloud2)
# install webshot
library(webshot)
webshot::install_phantomjs()

# Make the graph
my_graph <- wordcloud2(demoFreq, size = 1.5)

# save it in html
library("htmlwidgets")
saveWidget(my_graph, "./BLS_EP_files/tmp.html", selfcontained = F)

# and in png or pdf
webshot(
  "./BLS_EP_files/tmp.html",
  "./BLS_EP_files/fig_1.pdf",
  delay = 5,
  vwidth = 480,
  vheight = 480
)


################# Supports image files and graphic objects to be visualized
rm(list = ls())
require(magick)
require(ggplot2)
require(ggimage)
getwd()
dt <-
  read.csv("C:\\Users\\wane1\\Documents\\file\\sci\\aiep\\Annual_Production.csv")
dt <-
  read.csv("./BLS_EP_files/PubMed_Timeline_Results_by_Year_NM.csv")

dt$Year <- factor(dt$Year)
summary(dt)
windowsFonts(myFont = windowsFont("华文行楷"))
ggplot(data = dt, mapping = aes(x = factor(Year), y = Articles, group = 1)) +
  geom_line(colour = "#d5a478", linetype = 2, cex = 1.50) +
  geom_point(colour = "#d5a478") +
  geom_text(aes(label = Articles),position=position_dodge(width = 0.9),size = 4.5,vjust = -0.35)+
  xlab("发表年份") +
  ylab("发文量") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 50, 5), limits = c(0, 50)) +
  geom_bar(fill = "steelblue", stat = "identity", width = 0.5, position = position_dodge(0.6)) +
  theme_classic()+  theme(axis.text = element_text(size = 15,face = 'bold'))
ggsave()

library(showtext)
showtext_auto(enable = TRUE)
font_add(family ="YaHei",regular ='msyh.ttc')


p <- ggplot(dt, aes(x = Year, y = Articles, group = group)) +
  geom_line(aes(color = group)) +
  theme_classic() +
  ylab("Number of Publications Per Year") +
  geom_point(aes(color = group, shape = group)) +
  scale_color_discrete(
    name = " ",
    breaks = c("Radiology", "NM"),
    labels = c("AI within Radiology", "AI within Nuclear Medicine")
  ) +
  scale_shape_discrete(
    name = " ",
    breaks = c("Radiology", "NM"),
    labels = c("AI within Radiology", "AI within Nuclear Medicine")
  ) +
  theme(
    legend.justification = c(0.05, 1),
    legend.position = c(0.05, 1)
  )
# theme(legend.title=element_blank())
p
p + ggimage::geom_image()

img <- "./BLS_EP_files/pet_clinical-100_2.jpg"
ggbackground(p, img)
ggbackground(p, img, alpha = 0.00009) # color="steelblue"
ggbackground(
  p,
  img,
  image_fun = function(x) {
    image_negate(image_convolve(x, "DoG:0,10,10"))
  }
)
# Use custom color palettes
p + scale_color_manual(values = c("#E69F00", "#56B4E9"))
# Use brewer color palettes
p + scale_color_brewer(palette = "Dark2")
# Use grey scale
p + scale_color_grey() + theme_classic()
####### Add Background Image to ggplot2
library(ggpubr)
library(jpeg)
img <- readJPEG("./BLS_EP_files/SIGNA-PET-MR.jpg")
ggplot(dt, aes(x = Year, y = Publications, group = group)) +
  background_image(img) +
  geom_line(aes(color = group)) +
  theme_classic() +
  ylab("Number of Publications Per Year") +
  geom_point(aes(color = group, shape = group)) +
  scale_color_discrete(
    name = " ",
    breaks = c("Radiology", "NM"),
    labels = c("AI within Radiology", "AI within Nuclear Medicine")
  ) +
  scale_shape_discrete(
    name = " ",
    breaks = c("Radiology", "NM"),
    labels = c("AI within Radiology", "AI within Nuclear Medicine")
  ) +
  theme(
    legend.justification = c(0.05, 1),
    legend.position = c(0.05, 1)
  )
