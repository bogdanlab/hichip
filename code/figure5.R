# Figure 5 and analysis
# meta_analysis_functional_data.R
#######################################
# functions and libraries
library(metafor)
library(dplyr)
library(ggplot2) # plot
library(ggpubr) # plot

# https://stats.stackexchange.com/questions/30394/how-to-perform-two-sample-t-tests-in-r-by-inputting-sample-statistics-rather-tha
# equal.variance: whether or not to assume equal variance. Default is FALSE. 
t.test <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
    if( equal.variance==FALSE )
    {
        se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        # welch-satterthwaite df
        df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    } else
    {
        # pooled standard deviation, scaled by the sample sizes
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
        df <- n1+n2-2
    }
    t <- (m1-m2-m0)/se
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))
    names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
    return(dat)
}
#######################################
#######################################
# This is from Sandor's raw data "HiChip_KRAB_validation_raw_data.xlsx", Expression_AR for each of the Replicates 1,2,3 (T,U,V) columns. 
#######################################

ar1 = data.frame(OR = c(0.902, 1.031, 1.067, 0.649, 0.553, 0.509, 0.093, 0.025, 0.022, 0.033, 0.016, 0.009, 0.025, 0.011, 0.005), studlab = c(rep("CTRL", 3), rep("gRNA1",3),  rep("gRNA2",3), rep("gRNA3",3), rep("gRNA4",3)), tec_rep = rep(c("Rep1", "Rep2", "Rep3"), 5), biol_rep = "Rep1")
ar2 = data.frame(OR = c(0.99, 0.95, 1.06, 0.924, 0.936, 0.919, 0.701, 0.673, 0.678, 0.706, 0.743, 0.762, 0.798, 0.745, 0.711), studlab = c(rep("CTRL", 3), rep("gRNA1",3),     rep("gRNA2",3), rep("gRNA3",3), rep("gRNA4",3)), tec_rep = rep(c("Rep1", "Rep2", "Rep3"), 5), biol_rep = "Rep2")
ar = rbind.data.frame(ar1, ar2)

ar$logOR = log(ar$OR)

data = ar %>%
  group_by(studlab, biol_rep) %>%
  summarize(yi=mean(logOR), vi = (sd(logOR))^2)  %>% data.frame


res = data.frame()
reps = names(table(data$studlab))
for (r in reps) {
x = data[data$studlab==r,]
m = rma(yi, vi, data=x, method="FE")
m_df = data.frame(study = r, TE.fixed = m$beta, seTE.fixed = m$se, upper.fixed= m$ci.ub, lower.fixed=m$ci.lb, zval.fixed = m$zval, pval.fixed=m$pval)
res = rbind.data.frame(res, m_df)
}

res$TE.fixed= exp(res$TE.fixed)
res$upper.fixed = exp(res$upper.fixed)
res$lower.fixed = exp(res$lower.fixed)


# In the data frame "res": These are the results from meta-analyzing across two biological samples each of the CNTRL, gRNA1, gRNA2, gRNA3, gRNA4 that we report in Supplementary Table S11
#######################################
# We run Student's t-test pairwise between CNTRL and gRNA with equal.variance
# Add to the results in Table S11
#######################################


df = data.frame()
groups = as.character(res$study)
groups = groups[!(groups %in% "CTRL")]
eff_ctrl = res$TE.fixed[res$study=="CTRL"]
se_ctrl = res$seTE.fixed[res$study=="CTRL"]

stat.test = data.frame()
for (i in groups) {
eff_i = res$TE.fixed[res$study==i]
se_i = res$seTE.fixed[res$study==i]

t= c(eff_ctrl, eff_i, (se_ctrl*sqrt(2)), (se_i*sqrt(2)), "CTRL", i)
print(t)
stat.test = rbind.data.frame(stat.test, t(t))

x = t.test2( eff_ctrl, eff_i, (se_ctrl*sqrt(2)), (se_i*sqrt(2)), 2, 2, equal.variance=T)
x = data.frame(t(x))
df = rbind.data.frame(df, x)
}
stat.test$p.value = df$p.value
stat.test$p.value = signif(stat.test$p.value,1)

names(stat.test)[5] = "group1"
names(stat.test)[6] = "group2"

new.row = df[1,]
new.row[]=NA
df = rbind.data.frame(new.row, df)


res_all = cbind.data.frame(res, df)

#######################################
# Plot
#######################################
# g = ggplot(res_all, aes(x = study, y=TE.fixed)) + geom_bar(position=position_dodge(),  width=0.4, stat="identity",colour="black", size=.3) +         geom_errorbar(aes(ymin=TE.fixed-seTE.fixed, ymax=TE.fixed+seTE.fixed),size=.3, width=.2, position=position_dodge(.9)) + xlab("") + ylab("") + theme_bw()


pad = 0.01
label.df <- data.frame(study = res_all$study,
                       TE.fixed = res_all$TE.fixed+res_all$seTE.fixed+pad, pval = res_all$p.value)

label.df$pval = signif(label.df$pval, 1)
label.df$sig = cut(label.df$pval,breaks = c(-0.1, 0.0001, 0.001, 0.01, 0.05, 1),labels = c("****", "***", "**", "*", ""))
label.df$sig = as.character(label.df$sig)
label.df$sig[is.na(label.df$sig)] = ""
g1 = g + geom_text(data = label.df, label = label.df$sig)

library(ggpubr)
my_comparisons = list( c("CTRL", "gRNA1"), c("CTRL", "gRNA2"), c("CTRL", "gRNA3"), c("CTRL", "gRNA4") )

stat.test$sig = label.df$sig[-1]
stat.test$p.value_sig = paste(stat.test$p.value, stat.test$sig, sep="")
p = g + stat_pvalue_manual(
    data = stat.test, label = "p.value_sig", #"p.value",
    y.position = c(1.1,1.2,1.3,1.4)
    )

p1 = p + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

