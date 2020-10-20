setwd("D:/Faculdade/pibit/PIBIT/lidiane")

library(devtools)
install_github("ujjwalkarn/xda")
library(xda)
library(ggplot2)
library(forcats)



#install and load required packages
# BP >> CP >> MF
BP <- c("response to organic substance (GO:0010033)", 
        "response to chemical (GO:0042221)", 
        "cellular response to chemical stimulus (GO:0070887)", 
        "response to stimulus (GO:0050896)", 
        "regulation of biological quality (GO:0065008)", 
        "cell communication (GO:0007154)", 
        "cellular response to stimulus (GO:0051716)", 
        "signaling (GO:0023052)", 
        "signal transduction (GO:0007165)", 
        "G protein-coupled receptor signaling pathway (GO:0007186)")

CP <- c("integral component of plasma membrane (GO:0005887)", 
        "plasma membrane part (GO:0044459)", 
        "plasma membrane (GO:0005886)", 
        "plasma membrane region (GO:0098590)", 
        "synapse part (GO:0044456)", 
        "membrane part (GO:0044425)", 
        "synaptic membrane (GO:0097060)", 
        "membrane (GO:0016020)", 
        "GABA-A receptor complex (GO:1902711)", 
        "integral component of membrane (GO:0016021)")

MF <- c("molecular transducer activity (GO:0060089)", 
        "signaling receptor activity (GO:0038023)", 
        "transmembrane signaling receptor activity (GO:0004888)", 
        "neurotransmitter receptor activity (GO:0030594)", 
        "catalytic activity (GO:0003824)", 
        "phosphoric diester hydrolase activity (GO:0008081)", 
        "3',5'-cyclic-nucleotide phosphodiesterase activity (GO:0004114)", 
        "3',5'-cyclic-AMP phosphodiesterase activity (GO:0004115)", 
        "G protein-coupled amine receptor activity (GO:0008227)", 
        "transmembrane receptor protein kinase activity (GO:0019199)")

bp <- c(66, 74, 61, 94, 68, 79, 85, 78, 74, 42)
cp <- c(44, 55, 65, 28, 23, 69, 15, 77, 6, 57)
mf <- c(48, 45, 40, 15, 68, 12, 8, 7, 9, 10)

df <- data.frame(Term = rep(c(BP, CP, MF),1),
                Count = c(bp, cp, mf),
                GO = rep.int(c("BP","CP","MF"),c(10,10,10)))


# tem que inverter a ordem da lista antes (grupo por grupo)
# Reorder following the value of another column:
df %>%
  mutate(Term = fct_relevel(Term, 
                            "transmembrane receptor protein kinase activity (GO:0019199)",
                            "G protein-coupled amine receptor activity (GO:0008227)", 
                            "3',5'-cyclic-AMP phosphodiesterase activity (GO:0004115)", 
                            "3',5'-cyclic-nucleotide phosphodiesterase activity (GO:0004114)", 
                            "phosphoric diester hydrolase activity (GO:0008081)", 
                            "catalytic activity (GO:0003824)", 
                            "neurotransmitter receptor activity (GO:0030594)", 
                            "transmembrane signaling receptor activity (GO:0004888)", 
                            "signaling receptor activity (GO:0038023)", 
                            "molecular transducer activity (GO:0060089)",
                            "integral component of membrane (GO:0016021)",
                            "GABA-A receptor complex (GO:1902711)", 
                            "membrane (GO:0016020)", 
                            "synaptic membrane (GO:0097060)", 
                            "membrane part (GO:0044425)", 
                            "synapse part (GO:0044456)", 
                            "plasma membrane region (GO:0098590)", 
                            "plasma membrane (GO:0005886)", 
                            "plasma membrane part (GO:0044459)", 
                            "integral component of plasma membrane (GO:0005887)",
                            "G protein-coupled receptor signaling pathway (GO:0007186)",
                            "signal transduction (GO:0007165)", 
                            "signaling (GO:0023052)", 
                            "cellular response to stimulus (GO:0051716)", 
                            "cell communication (GO:0007154)", 
                            "regulation of biological quality (GO:0065008)", 
                            "response to stimulus (GO:0050896)", 
                            "cellular response to chemical stimulus (GO:0070887)", 
                            "response to chemical (GO:0042221)", 
                            "response to organic substance (GO:0010033)")) %>%
  ggplot(aes(fill=GO, x=Term, y=Count)) + 
  geom_bar(stat = "identity", alpha=.8, width=.6) +
  scale_fill_brewer(palette = "Set2") +
  coord_flip() +
  xlab("") +
  theme_bw()
