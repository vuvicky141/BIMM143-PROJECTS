#' ---
#' title: "Class 05: Data Visualization and Graphs in R"
#' author: "Vicky Vu "
#' date: "January 23, 2020"
#' ---

# Class 5 : Data Visualization and Exploration in R

#  Customizing Plots 
## Creating A Project and getting data to plot 

#First download the zip file of data to work with. The file should be moved into the folder that you are working on.


plot(1:10, col="blue", typ="o")
# 2A line plot 
# Need to import/read input data fire first 
baby <- read.table("bimm143_05_rstats/weight_chart.txt", header=TRUE)

# Ploting data with line graph 
# A basic age v. weight scatterplot
plot(baby$Age, baby$Weight)

# A scatterplot to  lineplot type="o", "pch" changes the shape of plot 
# "cex" changes the size of the 
plot(baby$Age, baby$Weight, type="o", col="blue", pch=15, lwd=2)

# ylim=c(2.10) to change the range of the y axis 
plot(baby$Age, baby$Weight, type="o", col="blue", pch=15, lwd=2, ylim=c(2,10))

# To change the x-axis label xlab=
plot(baby$Age, baby$Weight, type="o", col="blue", cex=1, pch=15, lwd=2, ylim=c(2,10), xlab="Age", ylab="Weight (kg)", main="Age V. Weight")

# 2B Barplot
mouse <- read.table("bimm143_05_rstats/feature_counts.txt", header=TRUE, sep="\t")

# Plot barplot 
par(mar=c(5,11,2,1))
barplot(mouse$Count, horiz=TRUE, col="pink", xlab="Feature", names.arg=mouse$Feature,las=1)

#another example 
par(mar=c(1,4,5,2))
plot(1:10)

# 2C Histogram 


# Section 3: using color in plots 
# 3A: Providing color vectors 
#counts <- read.delim("bimm143_05_rstats/male_female_counts.txt") also works as a shortcut 
counts <- read.table("bimm143_05_rstats/male_female_counts.txt", header=TRUE, sep="\t")
barplot(counts$Count, names.arg = counts$Sample, col=rainbow(nrow(counts)))

# another plot of the same thing with differnt colors
barplot(counts$Count, names.arg = counts$Sample, col=c("red","blue","dark green"), las=1)
        
# 3B Coloring by Value 
#importing up_down_expression.txt
up_down <- read.delim("bimm143_05_rstats/up_down_expression.txt")

        










