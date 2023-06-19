# Transdiagnostic Cognitive Biases Meta-Analysis 
# Lavigne KM, Deng J, & Sauvé G
# 2023

# Install packages ----
packages <- c( "dplyr", "meta", "metafor", "netmeta", "readxl", "rgl", "reshape2", "writexl")

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

ipak(packages)

### Install dmetar
# if (!require("remotes")) {
#   install.packages("remotes")
# }
# remotes::install_github("MathiasHarrer/dmetar")
library(dmetar)

# Data Preparation ----

## Get data
setwd("C:/Users/katie/OneDrive - McGill University/projects/transdx-cog-bias-review/meta-analysis")
df <- read_excel("DataExtractionWide.xlsx")

## Filter & make numeric
df <- df %>% 
  filter(Include == "YES") %>%
  mutate_at(c("n1", "m1", "sd1", "perc1",
              "n2", "m2", "sd2", "perc2",
              "n3", "m3", "sd3", "perc3",
              "n4", "m4", "sd4", "perc4",
              "n5", "m5", "sd5", "perc5",
              "n6", "m6", "sd6", "perc6",
              "nCON", "mCON", "sdCON", "percCON",
              "n_cor", "corr"), as.numeric)
df$JBI_category <- ifelse(df$JBI_score >= .8, "L", "H")

## Generate pairwise comparisons
df$Citation_unique <- make.names(df$Citation, unique=TRUE)
pw <- pairwise(treat = list(g1, g2, g3, g4, g5, g6, gCON), 
               n = list(n1, n2, n3, n4, n5, n6, nCON),
               mean = list(m1, m2, m3, m4, m5, m6, mCON),
               sd = list(sd1, sd2, sd3, sd4, sd5, sd6, sdCON),
               studlab = Citation_unique, data = df, sm = "SMD")

## Reverse code
### Means
pw_rev <- pw %>% 
  mutate(mean1 = case_when(
    reverse_code == "YES" ~ -mean1,
    reverse_code == "NO" ~ mean1,
    TRUE ~ 0)) %>%
  mutate(mean2 = case_when(
    reverse_code == "YES" ~ -mean2,
    reverse_code == "NO" ~ mean2,
    TRUE ~ 0))

### Correlations
df_rev_corr <- df %>%
  mutate(corr = case_when(
    reverse_code == "YES" ~ -corr,
    reverse_code == "NO" ~ corr,
    TRUE ~ 0))
pw$esid <- ave(pw$Citation, pw$Citation, FUN=seq_along)

## Save outputs
### Full data
save(pw_rev, file="full-pairwise-data.R")
write_xlsx(pw_rev, "full-pairwise-data.xlsx")

### Control comparisons
con <- pw_rev %>%
  filter(treat1 != "CONLR") %>%
  filter(treat2 == "CON" | treat2 == "CONLR")
con <- con[!is.na(con$mean1), ]
con <- con[!is.na(con$mean2), ]
#con$esid <- ave(con$Citation, con$Citation, FUN=seq_along) # estimate ID for aggregate future testing
save(con, file="control-pairwise-data.R")

### Symptom associations
sx <- df_rev_corr[!is.na(df_rev_corr$corr), ]
save(sx, file="symptoms-data.R")

# M1: Traditional meta-analysis on all patients versus controls on all biases ---- 

## Overall meta-analysis w/ risk of bias subgroups # more appropriate to do multivariate here
m1 <- metacont(n1, mean1, sd1,
               n2, mean2, sd2, 
               comb.random = T,
               data = con, 
               sm = "SMD",
               studylab='Citation_unique',
               random = TRUE,
               subgroup = JBI_category, 
               tau.common = FALSE,
               title = "Transdiagnostic Cognitive Biases - Patient vs. control")

sink(file="M1-pt-vs-con.txt", type = "output")
m1
sink()

## Subgroup analyses
### Bias categories (ATT, MEM, INT)
sink(file="M1-pt-vs-con-bias-category.txt", type = "output")
update.meta(m1, 
            subgroup = bias_category, 
            tau.common = FALSE)
sink()

## Symptom-specificity of content (Symptom-Specific, Neutral, Mixed)
sink(file="M1-pt-vs-con-symptom-specificity.txt", type = "output")
update.meta(m1, 
            subgroup = content, 
            tau.common = FALSE)
sink()

## Funnel Plot
pdf(file="funnel-plot-meta.pdf", width=50, height=50)
funnel(m1)
dev.off()

# M2: Network meta-analysis on diagnostic groups ----

## merge ANXa & ANXb to AD / CONLR to CON / SZd SZnd SZh SZnh PSY to SZ / 
dat <- pw_rev %>%
  mutate_at(vars(c("treat1", "treat2")), ~ str_replace(., "ANXa","AD")) %>%
  mutate_at(vars(c("treat1", "treat2")), ~ str_replace(., "ANXb","AD")) %>%
  mutate_at(vars(c("treat1", "treat2")), ~ str_replace(., "CONLR","CON")) %>%
  mutate_at(vars(c("treat1", "treat2")), ~ str_replace(., "SZd","SZ")) %>%
  mutate_at(vars(c("treat1", "treat2")), ~ str_replace(., "SZnd","SZ")) %>%
  mutate_at(vars(c("treat1", "treat2")), ~ str_replace(., "SZh","SZ")) %>%
  mutate_at(vars(c("treat1", "treat2")), ~ str_replace(., "SZnh","SZ")) %>%
  mutate_at(vars(c("treat1", "treat2")), ~ str_replace(., "PSY","SZ"))

## Filter (to do - loop)
nma <- dat %>%
  filter(bias_category == "ATT")
nma$esid <- ave(nma$Citation, nma$Citation, FUN=seq_along)
nma$study <- paste(nma$Citation, nma$esid)
nma <- nma[nma$treat1 != nma$treat2, ]

### NMA
m.netmeta <- netmeta(TE = TE,
                     seTE = seTE,
                     treat1 = treat1,
                     treat2 = treat2,
                     studlab = study,
                     data = nma,
                     sm = "SMD",
                     fixed = TRUE,
                     random = TRUE,
                     reference.group = "CON",
                     details.chkmultiarm = TRUE,
                     sep.trts = " vs ")

sink(file=paste("M2-nma-summary-ATT.txt", sep=""), type = "output")
summary(m.netmeta)
sink()

# decomp.design(m.netmeta) # Used to double-check if a random effects model should be used - it should.

## Graph figure
  ### 2d
  netgraph(m.netmeta,
           cex=1.25,
           offset=0.025,
           labels = m.netmeta$trts)
  ### 3d
  #netgraph(m.netmeta, dim = "3d")

## Direct & indirect evidence
d.evidence <- direct.evidence.plot(m.netmeta)
plot(d.evidence)

## Effect estimate table
result.matrix <- m.netmeta$TE.random
#result.matrix <- round(result.matrix, 2)
#result.matrix[lower.tri(result.matrix, diag = FALSE)] <- NA
r <- result.matrix

# Heatmap 
# Dummy data
x <- paste(c("AD","AN","BDD","BN","CON","DID","GAD","HA","MDD","OCD","PD","SAD")) # ATT
x <- paste(c("AD","CON","HA","MDD","OCD","PD","SAD","SZ")) # MEM
x <- paste(c("AD","AN","BN","BP","CON","GAD","HA","MDD","OCD","PD","PTSD","SAD","SD","SZ")) # INT
y <- x

melted_cormat <- melt(r)

ggplot(melted_cormat, aes(Var1, ordered(Var2, levels = rev(sort(unique(Var2)))))) +   
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "darkred") +
  theme(text = element_text(size = 16)) +
  xlab("") +
  ylab("")

## Produce effect table
netleague <- netleague(m.netmeta, 
                       bracket = "(", # use round brackets
                       digits=2)      # round to two digits

## Save results (here: the ones of the fixed-effect model)
write.csv(netleague$fixed, paste("M2-nma-netgraph-", bias, "-netleague.csv", sep=""))

## Ranking
netrank(m.netmeta, common=FALSE, small.values = "bad")

## Forest
forest(m.netmeta, 
       reference.group = "CON",
       xlim = c(-1.5, 1.5),
       smlab = bias,
       drop.reference.group = TRUE,
       label.left = "Shows no bias",
       label.right = "Shows bias",
       labels = m.netmeta$trts)

## Heatmap
netheat(m.netmeta, random = TRUE)

## Net splitting
netsplit(m.netmeta)

## Network splitting - forest
netsplit(m.netmeta) %>% forest()

# M3: Traditional meta-analysis on symptom associations with cognitive biases ----
## Overall Symptoms
m3a <- metacor(corr, n_cor, 
              comb.random = T,
              data = sx, 
              sm = "ZCOR",
              studylab='Citation_unique',
              random = TRUE,
              subgroup = JBI_category, 
              tau.common = FALSE,
              title = "Transdiagnostic Cognitive Biases - Symptom Associations")
sink(file="M3-sx-assocs.txt", type = "output")
m3a
sink()

## Subgroup analyses
### Bias categories (ATT, MEM, INT)
sink(file="M3-sx-assocs-bias-category.txt", type = "output")
update.meta(m3a, 
            subgroup = bias_category, 
            tau.common = FALSE)
sink()

## Symptom-specificity of content (Symptom-Specific, Neutral, Mixed)
sink(file="M3-sx-assocs-symptom-specificity.txt", type = "output")
update.meta(m3a, 
            subgroup = content, 
            tau.common = FALSE)
sink()

## Symptom categories
sx_cat <- sx[!is.na(sx$hitop_superspectra), ]
m3b <- metacor(corr, n_cor, 
              comb.random = T,
              data = sx_cat, 
              sm = "ZCOR",
              studylab='Citation_unique',
              random = TRUE, 
              subgroup = JBI_category, 
              tau.common = FALSE,
              title = "Transdiagnostic Cognitive Biases - Symptom Associations")
sink(file="M3-sx-assocs-sx-categories.txt", type = "output")
m3b
sink()

## Symptom-specificity of content (Symptom-Specific, Neutral, Mixed)
sink(file="M3-sx-assocs-sx-categories-hitop.txt", type = "output")
update.meta(m3b, 
            subgroup = hitop_superspectra, 
            tau.common = FALSE)
sink()