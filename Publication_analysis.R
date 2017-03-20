###########
#Libraries#
###########
library(sqldf)

###########################
#Global settings and files#
###########################
options(stringsAsFactors=F)
wd <- c('/usr3/graduate/akram/metabolic_modeling/metabolic_modeling')
od <- c('./modeling_output_files/')
setwd(wd)

annotation_reactions <- read.delim(file="sybil_Recon21x_model_react.tsv") 

###########
#Functions#
###########
source('~/metabolic_modeling/metabolic_modeling/metabolic_modeling_functions_linux.R')

#####################################################
#Comparing brown vs white predictions using seahorse#
#####################################################
mat_ws <- read.csv(file=paste0(od, "model_seahorse_white_predictions.csv"), row.names=1)
mat_bs <- read.csv(file=paste0(od, "model_seahorse_brown_predictions.csv"), row.names=1)

react_names <- rownames(mat_ws)
sample_test <- vector(mode="numeric",length=nrow(mat_ws))
sample_means <- vector(mode="numeric",length=nrow(mat_ws))
mat_1 <- as.matrix(mat_bs)
mat_2 <- as.matrix(mat_ws)
bat_mean <- vector(mode="numeric",length=nrow(mat_ws))
wat_mean <- vector(mode="numeric",length=nrow(mat_ws))

t.test_v2 <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
} 

for(i in 1:length(sample_test)){
  sample_test[i] <- t.test_v2(abs(mat_1[i,]),abs(mat_2[i,]))
  sample_means[i] <- (abs(mean(mat_1[i,])) - abs(mean(mat_2[i,])))
  bat_mean[i] <- mean(mat_1[i,])
  wat_mean[i] <- mean(mat_2[i,])
}

sample_test[sample_test == "NaN"] <-1
sample_fdr <- p.adjust(sample_test, method="BY")
results_mat <- cbind(sample_test, sample_fdr, sample_means, bat_mean, wat_mean)
rownames(results_mat) <- react_names
colnames(results_mat) <- c("p", "FDR", "Difference", "BAT", "WAT")
results_mat <- results_mat[order(abs(results_mat[,"FDR"]), decreasing=F),]
results_mat <- bind_annotation(results_mat)
rownames(results_mat) <- results_mat[,"react_ids"]

###################################################################
#Comparing brown vs white predictions using seahorse and nova data#
###################################################################
mat_ws <- read.csv(file=paste0(od, "model_seahorse_nova_white_predictions.csv"), row.names=1)
mat_bs <- read.csv(file=paste0(od, "model_seahorse_nova_brown_predictions.csv"), row.names=1)

react_names <- rownames(mat_ws)
sample_test <- vector(mode="numeric",length=nrow(mat_ws))
sample_means <- vector(mode="numeric",length=nrow(mat_ws))
mat_1 <- as.matrix(mat_bs)
mat_2 <- as.matrix(mat_ws)
bat_mean <- vector(mode="numeric",length=nrow(mat_ws))
wat_mean <- vector(mode="numeric",length=nrow(mat_ws))

t.test_v2 <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
} 

for(i in 1:length(sample_test)){
  sample_test[i] <- t.test_v2(abs(mat_1[i,]),abs(mat_2[i,]))
  sample_means[i] <- (abs(mean(mat_1[i,])) - abs(mean(mat_2[i,])))
  bat_mean[i] <- mean(mat_1[i,])
  wat_mean[i] <- mean(mat_2[i,])
}

sample_test[sample_test == "NaN"] <-1
sample_fdr <- p.adjust(sample_test, method="BY")
results_mat <- cbind(sample_test, sample_fdr, sample_means, bat_mean, wat_mean)
rownames(results_mat) <- react_names
colnames(results_mat) <- c("p", "FDR", "Difference", "BAT", "WAT")
results_mat <- results_mat[order(abs(results_mat[,"FDR"]), decreasing=F),]
results_mat <- bind_annotation(results_mat)
rownames(results_mat) <- results_mat[,"react_ids"]

write.csv(results_mat, file=paste0(od, "model_predicted_differences_results_seahorse_nova.csv"))

###########################################################################
#Comparing brown vs white predictions using seahorse and nova data for FVA#
###########################################################################
load("/projectnb/bu-joslin/metabolic_modeling/Tseng_adipocyte_metabolic_modeling_revised_model_v11_wa.RData")
mat_fva_lb_w <- read.csv(file=paste0(od, "model_seahorse_nova_white_FVA_predictions_lb_matrix.csv"), row.names=1)
mat_fva_ub_w <- read.csv(file=paste0(od, "model_seahorse_nova_white_FVA_predictions_ub_matrix.csv"), row.names=1)
mat_fva_lb_b <- read.csv(file=paste0(od, "model_seahorse_nova_brown_FVA_predictions_lb_matrix.csv"), row.names=1)
mat_fva_ub_b <- read.csv(file=paste0(od, "model_seahorse_nova_brown_FVA_predictions_ub_matrix.csv"), row.names=1)

#Look for mutual exclusive regions
me_1 <- rowSums((mat_fva_lb_b + 1000) > (mat_fva_ub_w + 1000))
me_2 <- rowSums((mat_fva_lb_w + 1000) > (mat_fva_ub_b + 1000))
crit <- me_1 + me_2 >= 1
fva_mat <- cbind(mat_fva_lb_b, mat_fva_ub_b, mat_fva_lb_w, mat_fva_ub_w)
colnames(fva_mat) <- c("BAT_lb", "BAT_ub", "WAT_lb","WAT_ub")
fva_mat <- bind_annotation(fva_mat)
View(fva_mat)

#Remove the numerical artifacts
fva_mat <- fva_mat[rowSums(abs(fva_mat[,2:5])) > 1e-6,]
write.csv(fva_mat, file=paste0(od, "FVA_mat_nonoverlapping.csv"))

#For each mutually exclusive region, make a plot
# library(stringr)
# 
# for(i in 1:(nrow(fva_mat)-3)){
#   m2p <- fva_mat[i,1:5]
#   rec <- fva_mat[i,"react_ids"]
#   des <- fva_mat[i,"name"]
#   m2p <- melt(m2p)
#   m2p[,"variable"] <- c("BA", "BA", "WA", "WA")
#   m2p <- cbind(m2p, c("lower", "upper", "lower", "upper"))
#   colnames(m2p) <- c("react_ids", "Tissue", "value", "range")
#   m2p <- dcast(m2p, Tissue ~ range)
#   
#   ggplot(m2p, aes(x=Tissue))+ 
#     geom_linerange(aes(ymin=lower, ymax=upper, col=Tissue)) +
#     geom_point(aes(y=lower, col=Tissue)) +
#     geom_point(aes(y=upper, col=Tissue)) +
#     labs(x="", y="Flux [mmol/gDw/hr]", title=paste(rec, ":\n", str_wrap(des, width=30))) +
#     scale_colour_manual(name="Tissue", values=c("black", "darkgrey")) +
#     theme_bw() + 
#     theme(plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#           legend.position="none", plot.title= element_text(size=rel(.5)))+
#     theme(panel.border= element_blank())+
#     theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))
#   ggsave(filename=paste0(od2,"/output_files_Tv11/", rec,".png"), width=2, height=3)
# }
