#!/usr/bin/env Rscript 

library("tidyverse")
library("edgeR")

counts <- read.table("./counts/feature_counts", sep = "\t",
skip = 1, header = TRUE, stringsAsFactors = FALSE)

counts_table <- counts %>% 
	column_to_rownames("Geneid") %>%
	select(-Chr, -Start, -End, -Strand, -Length) %>%
	as.matrix

group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=counts_table, group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

#F-test
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)


