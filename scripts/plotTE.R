#! /path/to/Rscript
##############################
library(phytools)
library(geiger)
library(nlme)
library(evomap)
library(caper)
library(aplot)
library(ggplot2)
library(splitstackshape)
# "MiscFunctions.R" is from "Bishop, Peter J., Mark A. Wright, and Stephanie E. Pierce. "Whole-limb scaling of muscle mass and force-generating capacity in amniotes." PeerJ 9 (2021): e12574."
source("MiscFunctions.R")

##############################
# PSLG plot
tree<-read.nexus("phyldog.tree")
data<- read.csv("TE_Genomesize.tbl", sep = "\t")
data$Genome_Size = data$Genome_Size/1000000
data$TE_length = data$TE_length/1000000
data$Genespace_length = data$Genespace_length/1000000
check <- name.check(tree, data, data.names=data$Species)
check
comp.data<-comparative.data(tree, data, 
                            names.col="Species", vcv.dim=2, 
                            warn.dropped=TRUE)
comp.data
model1 <- pgls(TE_length~Genome_Size, lambda="ML", data=comp.data)
model2 <- pgls(Genespace_length~Genome_Size, lambda="ML", data=comp.data)
summary(model1)
summary(model2)

lambda_opt1 = unname(model1$param[2])
a1 = summary(model1)$coefficients[1,1]
b1 = summary(model1)$coefficients[2,1]
print(c(lambda_opt1, a1, b1))

lambda_opt2 = unname(model2$param[2])
a2 = summary(model2)$coefficients[1,1]
b2 = summary(model2)$coefficients[2,1]
print(c(lambda_opt2, a2, b2))

adj_r_squared1 <- summary(model1)$adj.r.squared[1,1]
adj_p_value1 <- summary(model1)$coefficients[2,4]
annotation1 <- data.frame(
  x = c(600),
  y = c(100),
  label1 = paste("Pearson's test corrected by the phylogery.", "\n", "adj.r2 = ", adj_r_squared1, "\n","p value = ", adj_p_value1,sep = "")
)

adj_r_squared2 <- summary(model2)$adj.r.squared[1,1]
adj_p_value2 <- summary(model2)$coefficients[2,4]
annotation2 <- data.frame(
  x = c(600),
  y = c(300),
  label2 = paste("Pearson's test corrected by the phylogery.", "\n", "adj.r2 = ", adj_r_squared2, "\n","p value = ", adj_p_value2,sep = "")
)

tree1 = rescale(tree,"lambda",lambda_opt1)
pGLS_ci1 = gls.ci.mod(data$TE_length,data$Genome_Size,vcv(tree1))
pGLS_pi1 = gls.pi.mod(data$TE_length,data$Genome_Size,vcv(tree1),1)
meanPPE1 = mean.PPE(data$Genome_Size,data$TE_length,a1,b1)
PanAmniote_regression1 = cbind(a1,b1,
                                b1 - qt(0.975,length(tree$tip.label)-2)*summary(model1)$coefficients[2,2],
                                b1 + qt(0.975,length(tree$tip.label)-2)*summary(model1)$coefficients[2,2],
                                meanPPE1,lambda_opt1)
colnames(PanAmniote_regression1) = c("intercept","slope","slope_lower_CI","slope_upper_CI","mean_PPE","lambda_opt")
query_PI1 = gls.predict(data$TE_length,data$Genome_Size,vcv(tree1),1,log10(Xq))

tree2 = rescale(tree,"lambda",lambda_opt2)
pGLS_ci2 = gls.ci.mod(data$Genespace_length,data$Genome_Size,vcv(tree2))
pGLS_pi2 = gls.pi.mod(data$Genespace_length,data$Genome_Size,vcv(tree2),1)
meanPPE2 = mean.PPE(data$Genome_Size,data$Genespace_length,a2,b2)
PanAmniote_regression2 = cbind(a2,b2,
                                b2 - qt(0.975,length(tree$tip.label)-2)*summary(model1)$coefficients[2,2],
                                b2 + qt(0.975,length(tree$tip.label)-2)*summary(model1)$coefficients[2,2],
                                meanPPE2,lambda_opt2)
colnames(PanAmniote_regression2) = c("intercept","slope","slope_lower_CI","slope_upper_CI","mean_PPE","lambda_opt")
query_PI2 = gls.predict(data$Gene_length,data$Genome_Size,vcv(tree2),1,log10(Xq))

plotData = data
plotLine1 = data.frame(x=c(min(data$Genome_Size),max(data$Genome_Size)),y=c((a1+b1*min(data$Genome_Size)),(a1+b1*max(data$Genome_Size))))
plotCI1 = data.frame(x=pGLS_ci1$CI.plot$X,lower=pGLS_ci1$CI.plot$Lower2.5,upper=pGLS_ci1$CI.plot$Upper2.5)
plotPI1 = data.frame(x=pGLS_pi1$PI.plot$X,lower=pGLS_pi1$PI.plot$Lower2.5,upper=pGLS_pi1$PI.plot$Upper2.5)

plotLine2 = data.frame(x=c(min(data$Genome_Size),max(data$Genome_Size)),y=c((a2+b2*min(data$Genome_Size)),(a2+b2*max(data$Genome_Size))))
plotCI2 = data.frame(x=pGLS_ci2$CI.plot$X,lower=pGLS_ci2$CI.plot$Lower2.5,upper=pGLS_ci2$CI.plot$Upper2.5)
plotPI2 = data.frame(x=pGLS_pi2$PI.plot$X,lower=pGLS_pi2$PI.plot$Lower2.5,upper=pGLS_pi2$PI.plot$Upper2.5)

p1 <- ggplot() + 
    geom_point(data=plotData, aes(x=Genome_Size,y=TE_length), color="black", size=2, shape=19) +
    geom_line(data=plotLine1, aes(x=x, y=y), color="blue", size=1) +
    geom_ribbon(data=plotCI1,aes(x=x, ymin=lower, ymax=upper), fill=rgb(0.5,0.5,0.5), alpha=0.6) +
    xlab("Genome_Size (Mb)") + ylab("TE_length (Mb)") +
    geom_text(data=annotation1, aes( x=x, y=y, label=label1)) +
    theme_classic()


p2 <- ggplot() + 
    geom_point(data=plotData, aes(x=Genome_Size,y=Genespace_length), color="black", size=2, shape=19) +
    geom_line(data=plotLine2, aes(x=x, y=y), color="blue", size=1) +
    geom_ribbon(data=plotCI2,aes(x=x, ymin=lower, ymax=upper), fill=rgb(0.5,0.5,0.5), alpha=0.6) +
    xlab("Genome_Size (Mb)") + ylab("Genespace_length (Mb)") +
    geom_text(data=annotation2, aes( x=x, y=y, label=label2)) +
    theme_classic()

Pglsplot1 <- plot_list(p2, p1, ncol=2, tag_levels = "A")
ggsave(file = "Pglsplot1.pdf", Pglsplot1, width = 200, height = 140, units = "mm")

tree<-read.nexus("phyldog.tree")
data<- read.csv("TE_Genomesize.tbl",sep="\t")
data$Genome_Size = data$Genome_Size/1000000
data$SINEs = data$SINEs/1000000
data$LTR_elements = data$LTR_elements/1000000
data$LINEs = data$LINEs/1000000
data$DNA = data$DNA/1000000

check <- name.check(tree, data, data.names=data$Species)
check
comp.data<-comparative.data(tree, data, 
                            names.col="Species", vcv.dim=2, 
                            warn.dropped=TRUE)
comp.data
model3 <- pgls(SINEs~Genome_Size, lambda="ML", data=comp.data)
model4 <- pgls(LTR_elements~Genome_Size, lambda="ML", data=comp.data)
model5 <- pgls(LINEs~Genome_Size, lambda="ML", data=comp.data)
model6 <- pgls(DNA~Genome_Size, lambda="ML", data=comp.data)
summary(model3)
summary(model4)
summary(model5)
summary(model6)

lambda_opt3 = unname(model3$param[2])
a3 = summary(model3)$coefficients[1,1]
b3 = summary(model3)$coefficients[2,1]
print(c(lambda_opt3, a3, b3))

lambda_opt4  = unname(model4$param[2])
a4 = summary(model4)$coefficients[1,1]
b4 = summary(model4)$coefficients[2,1]
print(c(lambda_opt4, a4, b4))

lambda_opt5 = unname(model5$param[2])
a5 = summary(model5)$coefficients[1,1]
b5 = summary(model5)$coefficients[2,1]
print(c(lambda_opt5, a5, b5))

lambda_opt6  = unname(model6$param[2])
a6 = summary(model6)$coefficients[1,1]
b6 = summary(model6)$coefficients[2,1]
print(c(lambda_opt6, a6, b6))

adj_r_squared3 <- summary(model3)$adj.r.squared[1,1]
adj_p_value3 <- summary(model3)$coefficients[2,4]
annotation3 <- data.frame(
  x = c(600),
  y = c(0),
  label3 = paste("adj.r2 = ", adj_r_squared3, "\n","p value = ", adj_p_value3,sep = "")
)

adj_r_squared4 <- summary(model4)$adj.r.squared[1,1]
adj_p_value4 <- summary(model4)$coefficients[2,4]
annotation4 <- data.frame(
  x = c(600),
  y = c(0),
  label4 = paste("adj.r2 = ", adj_r_squared4, "\n","p value = ", adj_p_value4,sep = "")
)

adj_r_squared5 <- summary(model5)$adj.r.squared[1,1]
adj_p_value5 <- summary(model5)$coefficients[2,4]
annotation5 <- data.frame(
  x = c(600),
  y = c(0),
  label5 = paste("adj.r2 = ", adj_r_squared5, "\n","p value = ", adj_p_value5,sep = "")
)

adj_r_squared6 <- summary(model6)$adj.r.squared[1,1]
adj_p_value6 <- summary(model6)$coefficients[2,4]
annotation6 <- data.frame(
  x = c(600),
  y = c(0),
  label6 = paste("adj.r2 = ", adj_r_squared6, "\n","p value = ", adj_p_value6,sep = "")
)

tree3 = rescale(tree,"lambda",lambda_opt3)
pGLS_ci3 = gls.ci.mod(data$SINEs,data$Genome_Size,vcv(tree3))
pGLS_pi3 = gls.pi.mod(data$SINEs,data$Genome_Size,vcv(tree3),1)
meanPPE3 = mean.PPE(data$Genome_Size,data$SINEs,a3,b3)
PanAmniote_regression3 = cbind(a3,b3,
                                b3 - qt(0.975,length(tree$tip.label)-2)*summary(model3)$coefficients[2,2],
                                b3 + qt(0.975,length(tree$tip.label)-2)*summary(model3)$coefficients[2,2],
                                meanPPE3,lambda_opt3)
colnames(PanAmniote_regression3) = c("intercept","slope","slope_lower_CI","slope_upper_CI","mean_PPE","lambda_opt")
query_PI3 = gls.predict(data$SINEs,data$Genome_Size,vcv(tree3),1,log10(Xq))

tree4 = rescale(tree,"lambda",lambda_opt4)
pGLS_ci4 = gls.ci.mod(data$LTR_elements,data$Genome_Size,vcv(tree4))
pGLS_pi4 = gls.pi.mod(data$LTR_elements,data$Genome_Size,vcv(tree4),1)
meanPPE4 = mean.PPE(data$LTR_elements,data$Gene_length,a4,b4)
PanAmniote_regression4 = cbind(a4,b4,
                                b4 - qt(0.975,length(tree$tip.label)-2)*summary(model4)$coefficients[2,2],
                                b4 + qt(0.975,length(tree$tip.label)-2)*summary(model4)$coefficients[2,2],
                                meanPPE4,lambda_opt4)
colnames(PanAmniote_regression4) = c("intercept","slope","slope_lower_CI","slope_upper_CI","mean_PPE","lambda_opt")
query_PI4 = gls.predict(data$LTR_elements,data$Genome_Size,vcv(tree4),1,log10(Xq))

tree5 = rescale(tree,"lambda",lambda_opt5)
pGLS_ci5 = gls.ci.mod(data$LINEs,data$Genome_Size,vcv(tree5))
pGLS_pi5 = gls.pi.mod(data$LINEs,data$Genome_Size,vcv(tree5),1)
meanPPE5 = mean.PPE(data$Genome_Size,data$LINEs,a5,b5)
PanAmniote_regression5 = cbind(a5,b5,
                                b5 - qt(0.975,length(tree$tip.label)-2)*summary(model5)$coefficients[2,2],
                                b5 + qt(0.975,length(tree$tip.label)-2)*summary(model5)$coefficients[2,2],
                                meanPPE5,lambda_opt5)
colnames(PanAmniote_regression5) = c("intercept","slope","slope_lower_CI","slope_upper_CI","mean_PPE","lambda_opt")
query_PI5 = gls.predict(data$LINEs,data$Genome_Size,vcv(tree5),1,log10(Xq))

tree6 = rescale(tree,"lambda",lambda_opt6)
pGLS_ci6 = gls.ci.mod(data$DNA,data$Genome_Size,vcv(tree6))
pGLS_pi6 = gls.pi.mod(data$DNA,data$Genome_Size,vcv(tree6),1)
meanPPE6 = mean.PPE(data$DNA,data$Gene_length,a6,b6)
PanAmniote_regression6 = cbind(a6,b6,
                                b6 - qt(0.975,length(tree$tip.label)-2)*summary(model6)$coefficients[2,2],
                                b6 + qt(0.975,length(tree$tip.label)-2)*summary(model6)$coefficients[2,2],
                                meanPPE6,lambda_opt6)
colnames(PanAmniote_regression6) = c("intercept","slope","slope_lower_CI","slope_upper_CI","mean_PPE","lambda_opt")
query_PI6 = gls.predict(data$DNA,data$Genome_Size,vcv(tree6),1,log10(Xq))

plotData = data
plotLine3 = data.frame(x=c(min(data$Genome_Size),max(data$Genome_Size)),y=c((a3+b3*min(data$Genome_Size)),(a3+b3*max(data$Genome_Size))))
plotCI3 = data.frame(x=pGLS_ci3$CI.plot$X,lower=pGLS_ci3$CI.plot$Lower2.5,upper=pGLS_ci3$CI.plot$Upper2.5)
plotPI3 = data.frame(x=pGLS_pi3$PI.plot$X,lower=pGLS_pi3$PI.plot$Lower2.5,upper=pGLS_pi3$PI.plot$Upper2.5)

plotLine4 = data.frame(x=c(min(data$Genome_Size),max(data$Genome_Size)),y=c((a4+b4*min(data$Genome_Size)),(a4+b4*max(data$Genome_Size))))
plotCI4 = data.frame(x=pGLS_ci4$CI.plot$X,lower=pGLS_ci4$CI.plot$Lower2.5,upper=pGLS_ci4$CI.plot$Upper2.5)
plotPI4 = data.frame(x=pGLS_pi4$PI.plot$X,lower=pGLS_pi4$PI.plot$Lower2.5,upper=pGLS_pi4$PI.plot$Upper2.5)

plotLine5 = data.frame(x=c(min(data$Genome_Size),max(data$Genome_Size)),y=c((a5+b5*min(data$Genome_Size)),(a5+b5*max(data$Genome_Size))))
plotCI5 = data.frame(x=pGLS_ci5$CI.plot$X,lower=pGLS_ci5$CI.plot$Lower2.5,upper=pGLS_ci5$CI.plot$Upper2.5)
plotPI5 = data.frame(x=pGLS_pi5$PI.plot$X,lower=pGLS_pi5$PI.plot$Lower2.5,upper=pGLS_pi5$PI.plot$Upper2.5)

plotLine6 = data.frame(x=c(min(data$Genome_Size),max(data$Genome_Size)),y=c((a6+b6*min(data$Genome_Size)),(a6+b6*max(data$Genome_Size))))
plotCI6 = data.frame(x=pGLS_ci6$CI.plot$X,lower=pGLS_ci6$CI.plot$Lower2.5,upper=pGLS_ci6$CI.plot$Upper2.5)
plotPI6 = data.frame(x=pGLS_pi6$PI.plot$X,lower=pGLS_pi6$PI.plot$Lower2.5,upper=pGLS_pi6$PI.plot$Upper2.5)
Pglsplot2 <- ggplot() + 
    geom_point(data=plotData, aes(x=Genome_Size,y=SINEs), color="#d7191c", size=2, shape=19) +
    geom_line(data=plotLine3, aes(x=x, y=y), color="#d7191c", size=1) +
    geom_ribbon(data=plotCI3,aes(x=x, ymin=lower, ymax=upper), fill="#d7191c", alpha=0.4) +
    xlab("Genome_Size (Mb)") + ylab("TE_length (Mb)") +
    geom_point(data=plotData, aes(x=Genome_Size,y=LTR_elements), color="#fdae61", size=2, shape=19) +
    geom_line(data=plotLine4, aes(x=x, y=y), color="#fdae61", size=1) +
    geom_ribbon(data=plotCI4,aes(x=x, ymin=lower, ymax=upper), fill="#fdae61", alpha=0.4) +
    xlab("Genome_Size (Mb)") + ylab("TE_length (Mb)") +
    geom_point(data=plotData, aes(x=Genome_Size,y=LINEs), color="#abdda4", size=2, shape=19) +
    geom_line(data=plotLine5, aes(x=x, y=y), color="#abdda4", size=1) +
    geom_ribbon(data=plotCI5,aes(x=x, ymin=lower, ymax=upper), fill="#abdda4", alpha=0.4) +
    xlab("Genome_Size (Mb)") + ylab("TE_length (Mb)") +
    geom_point(data=plotData, aes(x=Genome_Size,y=DNA), color="#2b83ba", size=2, shape=19) +
    geom_line(data=plotLine6, aes(x=x, y=y), color="#2b83ba", size=1) +
    geom_ribbon(data=plotCI6,aes(x=x, ymin=lower, ymax=upper), fill="#2b83ba", alpha=0.4) +
    xlab("Genome_Size (Mb)") + ylab("TE_length (Mb)") +
    theme_classic() +
    annotate("point", shape=19,x=300,y=46, color="#d7191c", size=2)+
    annotate("rect", fill="#d7191c",xmin = 280, xmax = 320, ymin = 45, ymax =47, alpha = .4)+
    annotate("text", x=510, y=46, label=paste("SINE",gsub("\n"," ", annotation3$label3)))+
    annotate("point", shape=19,x=300,y=62, color="#fdae61", size=2)+
    annotate("rect", fill="#fdae61",xmin = 280, xmax = 320, ymin = 61, ymax = 63, alpha = .4)+
    annotate("text", x=510, y=62, label=paste("LTR",gsub("\n"," ", annotation4$label4)))+
    annotate("point", shape=19,x=300,y=70, color="#abdda4", size=2)+
    annotate("rect", fill="#abdda4",xmin = 280, xmax = 320, ymin = 69, ymax =71, alpha = .4)+
    annotate("text", x=510, y=70, label=paste("LINE",gsub("\n"," ", annotation5$label5)))+
    annotate("point", shape=19,x=300,y=54, color="#2b83ba", size=2)+
    annotate("rect", fill="#2b83ba",xmin = 280, xmax = 320, ymin = 53, ymax =55, alpha = .4)+
    annotate("text", x=505, y=54, label=paste("DNA",gsub("\n"," ", annotation6$label6)))
ggsave(file = "Pglsplot2.pdf", Pglsplot2, width = 200, height = 140, units = "mm")

##############################
# Stack plot
colors <- c("#FDB462","#CCEBC5","#80B1D3","#F08080","#A1A1A1")
my_stack_plot_data <- read.csv("TE_Genomesize_stack.tbl", header=T, sep = "\t")
my_stack_plot_data$Species <- factor(my_stack_plot_data$Species, levels=c("Msep", "Mlor","Aigs", "Apla", "Harm", "Hvir",  "Sexi","Sfru", "Slit", "Tni"), ordered = TRUE)
Stackplot <- ggplot(my_stack_plot_data,aes(Species,length,fill=class),colour=colors) +
    geom_bar(stat="identity",position="stack") +
    ggtitle("bar_plot") +
    scale_fill_manual(values=colors) + 
    theme_bw() + 
    theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(colour="black"))
ggsave(file = "Stackplot.pdf", Stackplot, width = 200, height = 160, units = "mm")

##############################
# Bar plot   
colors<-c("#80B1D3", "#8DD3C7", "#FFFFB3" ,"#BEBADA" , "#FB8072")
my_bar_plot_data<-read.csv("Line_percent", header=T, sep = "\t")
my_bar_plot_data$Class3 = factor(my_bar_plot_data$Class3, levels=c('L1','R2','RTE','Jockey','I'))
LINEplot <- ggplot(my_bar_plot_data, aes(x = Species, y = per, fill = Class3)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values = colors)+
    theme_classic()
ggsave(file = "LINEplot.pdf", LINEplot, width = 200, height = 160, units = "mm")

##############################
# Violin plot
colors<-c("#FB8072", "#BEBADA", "#80B1D3","#8DD3C7","#FFFFB3")
my_violin_data<-read.table("kim2d",header=T)
my_violin_data$absLen <- ceiling(my_violin_data$absLen/10000000)
new_violin<-expandRows(my_violin_data, "absLen")
new_violin$Species <- factor(new_violin$Species, levels=c("Msep", "Mlor","Aigs", "Apla", "Harm", "Hvir",  "Sexi","Sfru", "Slit", "Tni"), ordered = TRUE)
ViolinPlot<-ggplot(new_violin, aes(x = Species, y = Divergence,fill = class3)) +
    geom_violin(trim=FALSE) + 
    scale_fill_manual(values=colors) +
    scale_y_log10(breaks = c(0,10, 20, 100)) + 
    facet_wrap (~class3) + 
    theme_bw() +
    # facet_grid (~class3) 
    theme(panel.grid=element_blank())
ggsave(file = "GGPLOT VIOLIN PLOT.pdf", ViolinPlot, width = 297, height = 110, units = "mm")
