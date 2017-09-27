#!/usr/bin/Rscript

rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load("ggplot2",
               "scales",
               "tidyr",
               "gplots",
               "RColorBrewer",
               "PoiClaClu",
               "org.Cannuum.eg.db",
               "AnnotationDbi",
               "grid",
               "gridExtra",
               "limma",
               "reshape2",
               "Hmisc",
               "tcltk",
               "R.utils",
               "NMF",
               "edgeR",
               "statmod",
               "data.table",
               "msa",
               "seqinr",
               "ape",
               "Gviz",
               "rtracklayer",
               "Rsamtools",
               "compiler")

enableJIT(3)

options(scipen=12)

fancy_exponents = function(input_label){

    output_label = c()

    for (ii in 1:length(input_label)) {
        l = input_label[[ii]]

        if(!(is.na(l))){

            if (l == 0) {
                number = 0

            } else {
                l_10 = log10(l)

                l_c = as.character(l_10)
                ten_exp = gsub("(?<=\\.).+$", "", l_c, perl=TRUE)
                ten_exp = gsub("\\.", "", ten_exp)

                if (l >= 1) {
                    dec_part = l_10 - as.numeric(ten_exp)
                    real_part = round(10 ^ dec_part, 2)

                } else if (l < 1) {
                    dec_part = l_10 - as.numeric(ten_exp)
                    real_part = round(10 ^ dec_part, 3)
                    real_part = real_num * 10

                    ten_exp = as.numeric(ten_exp) - 1

                }
                number = paste(as.character(real_part), "  %*%  ", "10^", as.character(ten_exp), sep="")
            }
            output_label[[ii]] = number
        } else {
            output_label[[ii]] = NA
        }
    }
    parse(text = output_label)
}

scientific_10 <- function(x) {
    replace = gsub("(e\\+)|(e)", " %*% 10^", scientific_format()(x))
    replace = gsub("0\\.0 %\\*% 10\\^00", "0", replace)
    parse(text=replace)
}

mult_format <- function() {
    function(x) format(100*x,digits = 2)
}

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


countToTpm <- function(counts, effLen){
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}

exactPoiCI <- function (X, conf.level=0.95) {
    alpha = 1 - conf.level
    upper <- 0.5 * qchisq(1-alpha/2, 2*X+2)
    lower <- 0.5 * qchisq(alpha/2, 2*X)
    return(c(lower, upper))
}

round_df <- function(x, digits) {
    #https://stackoverflow.com/a/29876220
    # round all numeric variables
    # x: data frame
    # digits: number of digits to round
    numeric_columns <- sapply(x, is.numeric)
    x[numeric_columns] <-  round(x[numeric_columns], digits)
    x
}

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)

    combined <- switch(position,
                       "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                              legend,
                                              ncol = 1,
                                              heights = unit.c(unit(1, "npc") - lheight, lheight)),
                       "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                             legend,
                                             ncol = 2,
                                             widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

    grid.newpage()
    grid.draw(combined)

    # return gtable invisibly
    invisible(combined)

}


eval_EBE = function(RVD_set, EBE, pwr_mat, relative_EBE="", first_x_perfect=4){

    RVD_list = unlist(strsplit(RVD_set, "-"))

    EBE_list = unlist(strsplit(EBE, ""))

    relative_EBE_list = unlist(strsplit(relative_EBE, ""))

    rel_power = 0
    ebe_power = 0
    strong_count = 0
    mismatch_count = 0
    first_x_perfect_failstatus = FALSE
    for(ii in seq(1:length(EBE_list))){
        nt = EBE_list[ii]
        rel_nt = relative_EBE_list[ii]
        rvd = RVD_list[ii]

        powers = pwr_mat[pwr_mat$RVD == rvd, ]

        good_matches = unlist(strsplit(as.character(powers$match), ""))

        if(nt %in% good_matches){
            strong_count = strong_count + powers$strength
        } else{
            mismatch_count = mismatch_count + 1

            if(ii <= first_x_perfect){
                first_x_perfect_failstatus = TRUE
            }
        }

        ebe_power = ebe_power + as.numeric(powers[, 2:5][nt])

        if(length(relative_EBE_list) == 0 | rel_nt == "N"){
            rel_power = rel_power + as.numeric(apply(powers[, 2:5], 1, FUN=max))
        } else{
            rel_power = rel_power + as.numeric(powers[, 2:5][rel_nt])
        }
    }

    final_power = (ebe_power / rel_power) * 100

    result = c(final_power, strong_count, mismatch_count, first_x_perfect_failstatus)

    return(result)
}


ticks = 1:10
ooms = 10^(-6:12)
breaks = as.vector(ticks %o% ooms)
show.labels <- c(F, F, F, F, F, F, F, F, F, T)
labels <- as.character(breaks * show.labels)
labels <- gsub("^0$", "phantom(0)", labels)
for(ii in 1:length(labels)){
    labels[[ii]] = gsub("^0.009001$", "10^-6", labels[[ii]])
    labels[[ii]] = gsub("^0.00001$", "10^-5", labels[[ii]])
    labels[[ii]] = gsub("^0.0001$", "10^-4", labels[[ii]])
    labels[[ii]] = gsub("^0.001$", "10^-3", labels[[ii]])
    labels[[ii]] = gsub("^0.01$", "10^-2", labels[[ii]])
    labels[[ii]] = gsub("^0.1$", "10^-1", labels[[ii]])
    labels[[ii]] = gsub("^1$", "10^0", labels[[ii]])
    labels[[ii]] = gsub("^10$", "10^1", labels[[ii]])
    labels[[ii]] = gsub("^100$", "10^2", labels[[ii]])
    labels[[ii]] = gsub("^1000$", "10^3", labels[[ii]])
    labels[[ii]] = gsub("^10000$", "10^4", labels[[ii]])
    labels[[ii]] = gsub("^100000$", "10^5", labels[[ii]])
    labels[[ii]] = gsub("^1000000$", "10^6", labels[[ii]])
    labels[[ii]] = gsub("^10000000$", "10^7", labels[[ii]])
    labels[[ii]] = gsub("^100000000$", "10^8", labels[[ii]])
    labels[[ii]] = gsub("^1000000000$", "10^9", labels[[ii]])
    labels[[ii]] = gsub("^10000000000$", "10^10", labels[[ii]])
    labels[[ii]] = gsub("^100000000000$", "10^11", labels[[ii]])
    labels[[ii]] = gsub("^1000000000000$", "10^12", labels[[ii]])
}


plot_theme = theme(plot.title=element_text(size=20, face="bold", margin=margin(10,0,10,0), hjust=0.5),
                   axis.text.x=element_text(size=16, color="black", margin=margin(7, 0, 0, 0)), axis.text.y=element_text(size=16, color="black", margin=margin(0, 7, 0, 0)),
                   axis.line.x=element_line(color="black", size=.75), axis.line.y=element_line(color="black", size=.75),
                   axis.title.x=element_text(size=18, face="bold", vjust=-0.4), axis.title.y=element_text(size=18, face="bold", vjust=0.4),
                   axis.ticks.length=unit(0.3, "cm"), axis.ticks.x=element_line(color="black"), axis.ticks.y=element_line(color="black"),
                   panel.background=element_blank(),
                   panel.grid.major=element_line(size=0.2, color="grey80"), panel.grid.minor=element_blank(),
                   strip.text.x=element_text(size=20, face="bold"), strip.text.y=element_text(size=20, face="bold"),
                   strip.background=element_rect(colour="white", fill="white"),
                   legend.title=element_text(size=18, face="bold"), legend.text=element_text(size=16), legend.background=element_rect(fill="white"),
                   legend.key=element_rect(fill="white", size=3), legend.key.size = unit(2, 'lines'))

#source("/home/sven/Documents/Uni/Masterthesis/Analyses/2017_ECW-123_Xcv_AvrBs3_AvrBs4/Bioconductor/eval_TSS.R")

now = now = format(Sys.time(), "%Y%m%d%H%M%S")

################################################################################
################################################################################
# DATA IMPORT edgeR
maindir = "/home/sven/Documents/Uni/Masterthesis/Analyses/2017_ECW-123_Xcv_AvrBs3_AvrBs4/Bioconductor/20170720015916"
genome_fasta = "/home/sven/Documents/Uni/Masterthesis/ECW-123_ncbi_CIDer/ECW-123_ncbi_CIDer_v7.fa"
genome_annot = "/home/sven/Documents/Uni/Masterthesis/ECW-123_ncbi_CIDer/ECW-123_ncbi_CIDer_annotation_v7.gff3"
#genome_annot = "/home/sven/Documents/Uni/Masterthesis/ECW-123_ncbi_CIDer/ECW-123_ncbi_CIDer_CLCrun1_annotation_v7_mod.gff3"
mod_annot = "/home/sven/Documents/Uni/Masterthesis/ECW-123_ncbi_CIDer/ECW-123_v7_mod.gff3"

rawdata_edge = read.table(file=paste(maindir, "/", "20170720015916_RfeatureCounts", sep=""), header=TRUE, sep="\t")

# create main output dirs
dir.create(paste(maindir, "/images", sep=""), showWarnings=FALSE)
dir.create(paste(maindir, "/mappings", sep=""), showWarnings=FALSE)
dir.create(paste(maindir, "/TALE_EBEs", sep=""), showWarnings=FALSE)
dir.create(paste(maindir, "/DEGs", sep=""), showWarnings=FALSE)

for(reg in c("up", "down")){
    dir.create(paste(maindir, "/", reg, sep=""), showWarnings=FALSE)
    dir.create(paste(maindir, "/", reg, "/TALE_EBEs", sep=""), showWarnings=FALSE)
    dir.create(paste(maindir, "/", reg, "/TALE_EBEs/results", sep=""), showWarnings=FALSE)
    dir.create(paste(maindir, "/", reg, "/TALE_EBEs/results/EBE_alignments", sep=""), showWarnings=FALSE)
    dir.create(paste(maindir, "/", reg, "/TALE_EBEs/results/EBE_sequences", sep=""), showWarnings=FALSE)
    dir.create(paste(maindir, "/", reg, "/TSS", sep=""), showWarnings=FALSE)
    dir.create(paste(maindir, "/", reg, "/TSS/images", sep=""), showWarnings=FALSE)
    dir.create(paste(maindir, "/", reg, "/TSS/statistics", sep=""), showWarnings=FALSE)
    dir.create(paste(maindir, "/", reg, "/DEG_illustrations_tables", sep=""), showWarnings=FALSE)
}

ggcol = ggplotColours(n=3)
lib_ggcol_o = c(rep(ggcol[1], 3), rep(ggcol[2], 3), rep(ggcol[3], 3))
lib_ggcol = c(rep(ggcol[1], 2), rep(ggcol[2], 3), rep(ggcol[3], 3))

gff_all = import.gff3(genome_annot)
#gff_all = gff_all[!(is.na(gff_all$gene)),]
# remove chloroplast and mitochondrial annotations
#gff_all = gff_all[-(seqnames(gff_all) == "NC_018552.1"), ]

if(genome_annot == "/home/sven/Documents/Uni/Masterthesis/ECW-123_ncbi_CIDer/ECW-123_ncbi_CIDer_CLC2_annotation_v7.gff3"){
    gff_all$gene = gff_all$Name
}

gff_cds = gff_all[gff_all$type == "CDS"]
gff_gene = gff_all[gff_all$type == "gene"]

chromosomes = as.vector(unique(as.character(seqnames(gff_all))))

names(rawdata_edge) = gsub("(.+CIDer\\.)|(\\.bam)", "", names(rawdata_edge))

#names(rawdata_edge) = gsub("trimmed.", "", names(rawdata_edge))
names(rawdata_edge) = gsub("_sorted", "", names(rawdata_edge))
names(rawdata_edge) = gsub("_24h", "", names(rawdata_edge))
names(rawdata_edge) = gsub("empty", "control", names(rawdata_edge))
names(rawdata_edge) = gsub("_a", "_A", names(rawdata_edge))
names(rawdata_edge) = gsub("_b", "_B", names(rawdata_edge))
names(rawdata_edge) = gsub("_c", "_C", names(rawdata_edge))

#rawdata_edge = rawdata_edge[-grep("MSTRG", rawdata_edge$X), ]

# Preparations for setting up the DEGList object later
data_names = factor(names(rawdata_edge)[3:ncol(rawdata_edge)], levels=names(rawdata_edge))

rawdata_edge = rawdata_edge[-c(grep("(TRNA)", rawdata_edge$GeneID)), ]
myids = as.matrix(rawdata_edge$GeneID)
#rawdata_edge$GeneID = gsub("(?<=[0-9])\\.(.|..)$", "", rawdata_edge$GeneID, perl=TRUE)

rawdata_edge_length = rawdata_edge

na_elements = is.na(rawdata_edge$GeneID)
rawdata_edge$GeneID[na_elements] = myids[na_elements]
rownames(rawdata_edge) <- NULL
rawdata_edge = rawdata_edge[!duplicated(rawdata_edge$GeneID), ]

row.names(rawdata_edge) = rawdata_edge$GeneID
rawdata_edge$GeneID = NULL

experiment_groups = factor(c("1_control", "1_control", "1_control", "2_AvrBs3", "2_AvrBs3", "2_AvrBs3" , "3_AvrBs4", "3_AvrBs4", "3_AvrBs4"))

rawdata_edge = cbind(rawdata_edge[,grep("control", names(rawdata_edge))], rawdata_edge[,grep("AvrBs3", names(rawdata_edge))], rawdata_edge[,grep("AvrBs4", names(rawdata_edge))])


# Generating density plot of unfiltered data
logcpms_unf = stack(as.data.frame(cpm(rawdata_edge, log=TRUE)))
names(logcpms_unf)[2] = "Samples"
dens_unf = ggplot(data=logcpms_unf, aes(x=values, col=Samples, group=Samples)) +
    ylab("Density [%]") +
    xlab(expression(bold(log[2]~CPM))) +
    scale_x_continuous(limits=c(-10, 15)) +
    scale_y_continuous(labels=mult_format()) +
    plot_theme +
    ggtitle("Unfiltered") +
    theme(legend.position='none') +
    stat_density(position="identity", geom="line", size=1.3)

keep = rowSums(cpm(rawdata_edge) >= 1) >= 2
prefilter_data = rawdata_edge[keep, ]

# Generating density plot of filtered data
logcpms_f = stack(as.data.frame(cpm(prefilter_data, log=TRUE)))
names(logcpms_f)[2] = "Samples"
dens_f = ggplot(data=logcpms_f, aes(x=values, col=Samples, group=Samples)) +
    ylab("") +
    xlab(expression(bold(log[2]~CPM))) +
    scale_x_continuous(limits=c(-10, 15)) +
    scale_y_continuous(labels=mult_format()) +
    plot_theme +
    ggtitle("Filtered") +
    stat_density(position="identity", geom="line", size=1.3)

pdf(paste(maindir, "/images/filtered_data.pdf", sep=""), 10, 5, pointsize=10, pagecentre=TRUE, colormodel="cmyk")
grid_arrange_shared_legend(dens_unf, dens_f, ncol=2, nrow=1,  position="right")
dev.off()

png(paste(maindir, "/images/filtered_data.png", sep=""), width=900, height=450, units="px", pointsize=30)
grid_arrange_shared_legend(dens_unf, dens_f, ncol=2, nrow=1,  position="right")
dev.off()


# Plotting library sizes
lib_names = c("control_A", "control_B", "control_C", "AvrBs3_A", "AvrBs3_B", "AvrBs3_C", "AvrBs4_A", "AvrBs4_B", "AvrBs4_C")
lib_reads = c(17572018, 18058932, 16002504, 16694139, 18865204, 16861970, 17824407, 16621312, 17554802)

lib_size_raw = data.frame(rep(lib_names, 2))
colnames(lib_size_raw) = c("sample")
lib_size_raw$number = c(colSums(rawdata_edge), lib_reads)
lib_size_raw$group = c(rep("counted", 9), rep("fragments", 9))

ggplot(data=NULL, aes(x=factor(sample, levels=lib_names), y=number)) +
    xlab("") +
    ylab("Library size [reads]") +
    scale_y_continuous(limits=c(0, 20000000), labels=scientific_10) +
    plot_theme +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none') +
    geom_bar(data=lib_size_raw[lib_size_raw$group=="fragments",], stat="identity",  fill="grey60") +
    geom_bar(data=lib_size_raw[lib_size_raw$group=="counted",], stat="identity",  fill=lib_ggcol_o)
setwd(paste(maindir, "/", "images/", sep=""))
ggsave("library_size.pdf", width=20, units="cm")
ggsave("library_size.png", width=20, units="cm", dpi=450)

# Poisson distance matrix
table_pois = prefilter_data
sample_poisdistances = PoissonDistance(t(table_pois))
sample_poisdistances_matrix = as.matrix(sample_poisdistances$dd)
rownames(sample_poisdistances_matrix) = colnames(table_pois)
colnames(sample_poisdistances_matrix) = colnames(table_pois)

# Poisson heatmap
datatype = c("pdf", "png")
for(type in datatype){
    if(type == "pdf"){
        pdf(paste(maindir, "/", "images/pois_samples.pdf", sep=""), 7, 5, pointsize=12, pagecentre=TRUE, colormodel="rgb")
    } else {
        png(paste(maindir, "/", "images/pois_samples.png", sep=""), width=1600, height=1000, units="px", pointsize=30)
    }

    heatmap.2(sample_poisdistances_matrix,
              main="", # heat map title
              ColSideColors=lib_ggcol_o,
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins=c(4.7,6.5),     # widens margins around plot
              col=colorRampPalette(rev(brewer.pal(9, "Blues")))(200),
              dendrogram="col",
              key.xlab="Poisson distance value",
              key.title=NA,
              adjRow=c(0.1,0.65),
              adjCol=c(0.9,0.4),
              srtCol=40,
              sepwidth=c(0.01,0.01),
              sepcolor="grey",
              colsep=1:ncol(sample_poisdistances_matrix),
              rowsep=1:nrow(sample_poisdistances_matrix),
              lhei=c(1.25, 4, 2),
              key.xtickfun = function() {
                  breaks = pretty(parent.frame()$breaks)
                  breaks = breaks[c(1,length(breaks))]
                  list(at = parent.frame()$scale01(breaks), labels = breaks)}
    )
    #par(lend = 1)
    #legend(0, 0.57,
    #       legend=c("control", "AvrBs3", "AvrBs4"),
    #       col=c(ggcol[1], ggcol[2], ggcol[3]),
    #       lty=1,
    #       lwd=10,
    #       cex=0.75)
    dev.off()
}
# http://stackoverflow.com/questions/23613866/heatmap-2-and-color-key-tick-marks
rm(table_pois)
rm(sample_poisdistances_matrix)
rm(sample_poisdistances)

# Generate MDS plot
mds = data.frame(plotMDS(prefilter_data[1:ncol(prefilter_data)], top=50000)$cmdscale.out)
mds$names = as.factor(gsub("^.+_", "", row.names(mds)))
ggplot(mds, aes(x=X1, y=X2, label=names, group=experiment_groups, shape=experiment_groups, colour=experiment_groups)) +
    xlab("Leading logFC dimension 1") +
    ylab("Leading logFC dimension 2") +
    plot_theme +
    theme(legend.position="right") +
    scale_shape_discrete(name="Conditions",
                         breaks=c("1_control", "2_AvrBs3", "3_AvrBs4"),
                         labels=c("control", "AvrBs3", "AvrBs4")) +
    scale_colour_discrete(name="Conditions",
                          breaks=c("1_control", "2_AvrBs3", "3_AvrBs4"),
                          labels=c("control", "AvrBs3", "AvrBs4")) +
    geom_text(aes(label=names), nudge_y=2, nudge_x=100, size=5) +
    geom_point(size=6) +
    theme(legend.position="top")
setwd(paste(maindir, "/", "images/", sep=""))
ggsave("MDS_plot.pdf", width=23, units="cm")
ggsave("MDS_plot.png", width=23, units="cm", dpi=450)

rm(prefilter_data)

rawdata_edge = rawdata_edge[, -(grep("control_C", names(rawdata_edge)))]
rawdata_edge_length = rawdata_edge_length[, -(grep("control_C", names(rawdata_edge_length)))]

experiment_groups = factor(c("1_control", "1_control", "2_AvrBs3", "2_AvrBs3", "2_AvrBs3", "3_AvrBs4", "3_AvrBs4", "3_AvrBs4"))
design = model.matrix(~experiment_groups, data=rawdata_edge)

# Setting up the edgeR DEGList object
edgeR_dge = DGEList(counts=rawdata_edge, group=experiment_groups, genes=rownames(rawdata_edge))
edgeR_dge$ENTREZ = mapIds(org.Cannuum.eg.db,
                          keys=as.vector(row.names(edgeR_dge)),
                          column="ENTREZID",
                          keytype="SYMBOL",
                          multiVals="first")
edgeR_dge$Full_Name = mapIds(org.Cannuum.eg.db,
                             keys=as.vector(row.names(edgeR_dge)),
                             column="GENENAME",
                             keytype="SYMBOL",
                             multiVals="first")
edgeR_dge$Lengths = rawdata_edge_length$Length

rm(rawdata_edge)
rawdata_edge_length$control_C = NA

# Filtering data
keep = rowSums(cpm(edgeR_dge) >= 1) >= 2
edgeR_dge = edgeR_dge[keep, , keep.lib.sizes=FALSE]

# Normalizing data
edgeR_dge = calcNormFactors(edgeR_dge, method="TMM")

# Generating plot of normalized data
logcpms = melt(as.data.frame(cpm(edgeR_dge, log=TRUE)))
ggplot(data=logcpms, aes(x=variable, y=value)) +
    xlab("") +
    ylab(expression(bold(log[2]~CPM))) +
    plot_theme +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none') +
    geom_boxplot(fill=lib_ggcol)
setwd(paste(maindir, "/", "images/", sep=""))
ggsave("normalized_counts.pdf", width=20, units="cm")
ggsave("normalized_counts.png", width=20, units="cm", dpi=450)

plotMD(edgeR_dge, column=1)
abline(h=0, col="red", lty=2, lwd=2)

# Estimate Dispersion
edgeR_dge = estimateDisp(edgeR_dge, design, robust=TRUE)

plotBCV(edgeR_dge)

# Create fit
fit = glmFit(edgeR_dge, design)

# Testing
lrt = glmLRT(fit, coef=2:3)

mean_cpm = NULL
mean_cpm_raw = as.data.frame(cpm(edgeR_dge, log=TRUE))
mean_cpm_raw$gene_symbols = row.names(mean_cpm_raw)
mean_cpm = data.frame(ID=row.names(mean_cpm_raw), Means=rowMeans(mean_cpm_raw[, grep("control", names(mean_cpm_raw))]))
mean_cpm = cbind(mean_cpm, data.frame(ID=row.names(mean_cpm_raw), Means=rowMeans(mean_cpm_raw[, grep("AvrBs3", names(mean_cpm_raw))])))
mean_cpm[,3] = NULL
mean_cpm = cbind(mean_cpm, data.frame(ID=row.names(mean_cpm_raw), Means=rowMeans(mean_cpm_raw[, grep("AvrBs4", names(mean_cpm_raw))])))
mean_cpm[,4] = NULL
names(mean_cpm) = c("gene_symbol", "mean_Control", "mean_AvrBs3", "mean_AvrBs4")

# Testing with LogFoldChange limit
t_2vs1 = glmTreat(fit, coef=2, lfc=1)
num_degs_2vs1 = summary(decideTestsDGE(t_2vs1))

t2vs1_df = as.data.frame(t_2vs1$table)
t2vs1_df = cbind(row.names(t2vs1_df), t2vs1_df)
t2vs1_df$unshrunk.logFC = NULL
names(t2vs1_df)[1] = "Gene_IDs"
t2vs1_df$treatment = "AvrBs3"
DEG = as.data.frame(decideTestsDGE(t_2vs1))
names(DEG)[1] = "DEG"
t2vs1_df = cbind(t2vs1_df, DEG)
t2vs1_df$logCPM = mean_cpm$mean_AvrBs3

t_3vs1 = glmTreat(fit, coef=3, lfc=1)
num_degs_3vs1 = summary(decideTestsDGE(t_3vs1))

t3vs1_df = as.data.frame(t_3vs1$table)
t3vs1_df = cbind(row.names(t3vs1_df), t3vs1_df)
t3vs1_df$unshrunk.logFC = NULL
names(t3vs1_df)[1] = "Gene_IDs"
t3vs1_df$treatment = "AvrBs4"
DEG = as.data.frame(decideTestsDGE(t_3vs1))
names(DEG)[1] = "DEG"
t3vs1_df = cbind(t3vs1_df, DEG)
t3vs1_df$logCPM = mean_cpm$mean_AvrBs4

t_both_df = rbind(t2vs1_df, t3vs1_df)


data_edge_length = merge(t3vs1_df, rawdata_edge_length, by.x="Gene_IDs", by.y="GeneID", all=FALSE)


t_2vs3 = glmTreat(fit, contrast=c(0, 1, -1), lfc=1)
num_degs_2vs3 = summary(decideTestsDGE(t_2vs3))

# p values
ggplot(data=t_both_df, aes(x=PValue)) +
    ylab("Count") +
    xlab("p values") +
    scale_y_continuous(breaks=seq(0, 15000, by=1000), labels=scientific_10) +
    scale_x_continuous(breaks=c(0, 0.5, 1),
                       labels=c("0", "0.5", "1")) +
    plot_theme +
    geom_histogram(col="black", fill="grey", binwidth=0.05, boundary=0) +
    facet_grid(. ~ treatment)
setwd(paste(maindir, "/", "images/", sep=""))
ggsave("p_vals.pdf", width=20, units="cm")
ggsave("p_vals.png", width=20, units="cm", dpi=450)


# log2 fold change
ggplot(data=t_both_df, aes(x=logFC)) +
    ylab("Count") +
    xlab(expression(bold(log[2]~fold~change))) +
    scale_y_log10(limits=c(1, NA), labels=parse(text=labels), breaks=breaks) +
    scale_x_continuous(limits=c(-12, 12),
                       breaks=c(-10, -5, 0, 5, 10),
                       labels=c("-10", "-5", "0", "5", "10")) +
    #annotation_logticks(base=10, sides="l") +
    plot_theme +
    geom_histogram(col="black", fill="grey", binwidth=1, center=0) +
    facet_grid(. ~ treatment)
setwd(paste(maindir, "/", "images/", sep=""))
ggsave("logFC.pdf", width=20, units="cm")
ggsave("logFC.png", width=20, units="cm", dpi=450)

data_edge_length = data_edge_length[, -c(2:6)]
data_edge_length = data_edge_length[, 1:(ncol(data_edge_length)-1)]
data_edge_TPM = data_edge_length
data_edge_TPM[, 3:ncol(data_edge_length)] = countToTpm(data_edge_length[, 3:ncol(data_edge_length)], data_edge_length$Length)

tpm = as.data.frame(countToTpm(data_edge_length[, 3:ncol(data_edge_length)], data_edge_length$Length))
tpm$Gene_Symbol = data_edge_length$Gene_IDs

# Coverage stuff
#
# gene_reads = GRanges(seqnames=as.character(seqnames(gene_info)), IRanges(window_start, window_end))
# what_info = c("rname", "pos")
# param = ScanBamParam(which=gene_reads, what=what_info, simpleCigar=TRUE)
#
# for(bam_name in bam_sets){
#     bamfile_path = paste(main_bams, "/", bam_name, sep="")
#
#     bamFile = BamFile(bamfile_path, index=paste(bamfile_path, ".", "bai", sep=""))
#
#     bamInfo = scanBam(bamFile, param=param)
#     reads_pos = as.numeric(unlist(bamInfo[[1]][2]))
# }


################################################################################
# Test result DF
test_res = cbind(t2vs1_df[1], as.data.frame(decideTestsDGE(t_2vs1)), as.data.frame(decideTestsDGE(t_3vs1)))
names(test_res) = c("Gene_Symbol", "AvrBs3", "AvrBs4")
test_res$Full_Name = mapIds(org.Cannuum.eg.db,
                            keys=as.vector(test_res$Gene_Symbol),
                            column="GENENAME",
                            keytype="SYMBOL",
                            multiVals="first")
test_res$AvrBs3_logFC = t2vs1_df$logFC
test_res$AvrBs4_logFC = t3vs1_df$logFC
test_res = merge(test_res, tpm, by="Gene_Symbol")
test_res$Control_meanTPM = rowMeans(test_res[, grep("control", names(test_res)[7:ncol(test_res)])+6])
test_res$AvrBs3_meanTPM = rowMeans(test_res[, grep("AvrBs3", names(test_res)[7:ncol(test_res)])+6])
test_res$AvrBs4_meanTPM = rowMeans(test_res[, grep("AvrBs4", names(test_res)[7:ncol(test_res)])+6])
test_res_full = cbind(test_res[1], test_res[4], test_res[2:3], test_res[5:6], test_res[15:17])
test_res_full = round_df(test_res_full, 3)

test_res_full_long_h = test_res_full
test_res_full_long_h = test_res_full_long_h[-grep("(Full_Name)|(Control)", names(test_res_full_long_h))]
names(test_res_full_long_h) = gsub("AvrBs3$", "AvrBs3_DEG", names(test_res_full_long_h))
names(test_res_full_long_h) = gsub("AvrBs4$", "AvrBs4_DEG", names(test_res_full_long_h))

test_res_full_long = data.frame(c(rep("AvrBs3", nrow(test_res_full_long_h)), rep("AvrBs4", nrow(test_res_full_long_h))))
test_res_full_long$DEG = c(as.matrix(test_res_full_long_h$AvrBs3_DEG), as.matrix(test_res_full_long_h$AvrBs4_DEG))
test_res_full_long$logFC = c(as.matrix(test_res_full_long_h$AvrBs3_logFC), as.matrix(test_res_full_long_h$AvrBs4_logFC))
test_res_full_long$meanTPM = c(as.matrix(test_res_full_long_h$AvrBs3_meanTPM), as.matrix(test_res_full_long_h$AvrBs4_meanTPM))
test_res_full_long$logTPM = log2(c(as.matrix(test_res_full_long_h$AvrBs3_meanTPM), as.matrix(test_res_full_long_h$AvrBs4_meanTPM)))
names(test_res_full_long) = c("TALE", "DEG", "logFC", "meanTPM", "logTPM")
test_res_full_long$logTPM[test_res_full_long$logTPM < -10] = -11
test_res_full_long$logTPM[test_res_full_long$logFC < -10] = -11

test_res_full_2 = test_res_full[,-2]
DEG_AvrBs3 = test_res_full_2[(test_res_full_2$AvrBs3 != 0),]
write.table(DEG_AvrBs3, file=paste(maindir, "/DEGs/degs_AvrBs3.tsx", sep=""), col.names=FALSE, row.names=FALSE, sep=" & ", quote=FALSE)

DEG_AvrBs4 = test_res_full_2[(test_res_full_2$AvrBs4 != 0),]
write.table(DEG_AvrBs4, file=paste(maindir, "/DEGs/degs_AvrBs4.tsx", sep=""), col.names=FALSE, row.names=FALSE, sep=" & ", quote=FALSE)

DEG_both = test_res_full_2[((test_res_full_2$AvrBs3 != 0) & (test_res_full_2$AvrBs4 != 0)), ]
write.table(DEG_both, file=paste(maindir, "/DEGs/degs_both.tsx", sep=""), col.names=FALSE, row.names=FALSE, sep=" & ", quote=FALSE)

DEG_both_oppositional = test_res_full_2[((test_res_full_2$AvrBs3 == 1) & (test_res_full_2$AvrBs4 == -1)) | ((test_res_full_2$AvrBs3 == -1) & (test_res_full_2$AvrBs4 == 1)),]
write.table(DEG_both_oppositional, file=paste(maindir, "/DEGs/degs_both_opp.tsx", sep=""), col.names=FALSE, row.names=FALSE, sep=" & ", quote=FALSE)



test_res_AvrBs3 = test_res[test_res$AvrBs3 == 1, ]
test_res_AvrBs3$meanControl = rowMeans(test_res_AvrBs3[, grep("control", names(test_res_AvrBs3))])
test_res_AvrBs3$meanAvrBs3 = rowMeans(test_res_AvrBs3[, grep("AvrBs3", names(test_res_AvrBs3))])
test_res_AvrBs3$meanAvrBs4 = rowMeans(test_res_AvrBs3[, grep("AvrBs4", names(test_res_AvrBs3))])

test_res_AvrBs4 = test_res[test_res$AvrBs4 == 1, ]
test_res_AvrBs4$meanControl = rowMeans(test_res_AvrBs4[, grep("control", names(test_res_AvrBs4))])
test_res_AvrBs4$meanAvrBs3 = rowMeans(test_res_AvrBs4[, grep("AvrBs3", names(test_res_AvrBs4))])
test_res_AvrBs4$meanAvrBs4 = rowMeans(test_res_AvrBs4[, grep("AvrBs4", names(test_res_AvrBs4))])

test_res_both = test_res[(test_res$AvrBs3 == 1) & (test_res$AvrBs4 == 1), ]
test_res_both$meanControl = rowMeans(test_res_both[, grep("control", names(test_res_both))])
test_res_both$meanAvrBs3 = rowMeans(test_res_both[, grep("AvrBs3", names(test_res_both))])
test_res_both$meanAvrBs4 = rowMeans(test_res_both[, grep("AvrBs4", names(test_res_both))])



write.table(test_res$Gene_Symbol, file=paste(maindir, "/", "TALE_EBEs/all_genes.txt", sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

upreg_all_genes = test_res[(test_res$AvrBs3 == 1) | (test_res$AvrBs4 == 1), ]
write.table(upreg_all_genes$Gene_Symbol, file=paste(maindir, "/", "TALE_EBEs/degs_all_up.txt", sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

downreg_all_genes = test_res[(test_res$AvrBs3 == -1) | (test_res$AvrBs4 == -1), ]
write.table(downreg_all_genes$Gene_Symbol, file=paste(maindir, "/", "TALE_EBEs/degs_all_down.txt", sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)


# MA plot
t_both_df_none = test_res_full_long[test_res_full_long$DEG == 0, ]
t_both_df_up = test_res_full_long[test_res_full_long$DEG == 1, ]
t_both_df_down = test_res_full_long[test_res_full_long$DEG == -1, ]
ggplot(data=NULL, aes(x=logTPM, y=logFC)) +
    ylab(expression(bold(log[2]~fold~change))) +
    xlab(expression(bold(log[2]~TPM))) +
    scale_y_continuous(breaks=seq(-14, 14, 2)) +
    scale_x_continuous(breaks=seq(-10, 14, 2)) +
    plot_theme +
    theme(plot.title=element_text(size=22, face="bold", margin=margin(10,0,10,0), hjust=0.5)) +
    geom_vline(xintercept=0, color="grey50") +
    geom_hline(yintercept=0, color="grey50") +
    geom_point(data=t_both_df_none, col="grey30") +
    geom_point(data=t_both_df_up, col="red", size=2) +
    geom_point(data=t_both_df_down, col="blue", size=2) +
    facet_grid(. ~ TALE)
setwd(paste(maindir, "/", "images/", sep=""))
ggsave("MA_plot.pdf", width=30, height=15, units="cm")
ggsave("MA_plot.png", width=25, units="cm", dpi=450)


# AvrBs3
upreg_AvrBs3_genes = test_res[test_res$AvrBs3 == 1, ]
#upreg_AvrBs3_genes = upreg_AvrBs3_genes[with(upreg_AvrBs3_genes, order(-AvrBs3_logFC)), ]
#upreg_AvrBs3_genes = upreg_AvrBs3_genes[c(1:25), ]

write.table(upreg_AvrBs3_genes$Gene_Symbol , file=paste(maindir, "/up/", "TALE_EBEs/degs_AvrBs3_up.txt", sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

downreg_AvrBs3_genes = test_res[test_res$AvrBs3 == -1, ]
write.table(downreg_AvrBs3_genes$Gene_Symbol , file=paste(maindir, "/down/", "TALE_EBEs/degs_AvrBs3_down.txt", sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

AvrBs3_small_logFC = test_res[, c(1, 5)]
AvrBs3_small_logFC$AvrBs3_logFC = 2^AvrBs3_small_logFC$AvrBs3_logFC
names(AvrBs3_small_logFC) = c("Gene_Symbol", "AvrBs3_logFC")


# AvrBs4
upreg_AvrBs4_genes = test_res[test_res$AvrBs4 == 1, ]
#upreg_AvrBs4_genes = upreg_AvrBs4_genes[with(upreg_AvrBs4_genes, order(-AvrBs4_logFC)), ]
#upreg_AvrBs4_genes = upreg_AvrBs4_genes[c(1:25), ]

write.table(upreg_AvrBs4_genes$Gene_Symbol , file=paste(maindir, "/up/", "TALE_EBEs/degs_AvrBs4_up.txt", sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

downreg_AvrBs4_genes = test_res[test_res$AvrBs4 == -1, ]
write.table(downreg_AvrBs4_genes$Gene_Symbol , file=paste(maindir, "/down/", "TALE_EBEs/degs_AvrBs4_down.txt", sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

AvrBs4_small_logFC = test_res[, c(1, 6)]
AvrBs4_small_logFC$AvrBs4_logFC = 2^AvrBs4_small_logFC$AvrBs4_logFC
names(AvrBs4_small_logFC) = c("Gene_Symbol", "AvrBs4_logFC")


# Both
reg_both_up = test_res[test_res$AvrBs3 == 1 & test_res$AvrBs4 == 1, ]
#reg_both_up = reg_both_up[with(reg_both_up, order(-AvrBs3_logFC)), ]
#reg_both_up = reg_both_up[c(1:25), ]

write.table(reg_both_up$Gene_Symbol , file=paste(maindir, "/up/", "TALE_EBEs/degs_both_up.txt", sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

reg_both_down = test_res[test_res$AvrBs3 == -1 & test_res$AvrBs4 == -1, ]
write.table(reg_both_up$Gene_Symbol , file=paste(maindir, "/down/", "TALE_EBEs/degs_both_down.txt", sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

# Venn Diagram
names(test_res)[2:3] = c("AvrBs3", "AvrBs4")
test_res_venn = test_res[2:3]
test_res_venn = abs(test_res_venn)

for(type in datatype){
    if(type == "pdf"){
        pdf(paste(maindir, "/images/venn.pdf", sep=""), 7, 5, pointsize=10, pagecentre=TRUE, colormodel="rgb")
    } else {
        png(paste(maindir, "/images/venn.png", sep=""), width=800, height=600, units="px", pointsize=16)
    }
    vennDiagram(test_res_venn, circle.col = c(ggcol[2], ggcol[3]))
    dev.off()
}

# Heatmap of DEGs
test_res_long = gather(test_res, key="treatment", value="DEG_flag", c(2, 3))
t_both_df = cbind(t_both_df, test_res_long[2])

test_res_long[3:6] = NULL
test_res_long$logFC = rbind(as.matrix(test_res$AvrBs3_logFC), as.matrix(test_res$AvrBs4_logFC))


t_both_df_deg = t_both_df[(t_both_df$DEG != 0), ]
deg_genes = as.vector(unique(t_both_df_deg$Gene_IDs))
deg_df = data.frame()
for(element in deg_genes){
    tmp = t_both_df[grep(paste("^", element, "$", sep=""), t_both_df$Gene_IDs), ]
    deg_df = rbind(deg_df, tmp)
}
AvrBs3_deg_df = deg_df[(deg_df$treatment == "AvrBs3"), ]
AvrBs4_deg_df = deg_df[(deg_df$treatment == "AvrBs4"), ]
degs_heatmap = as.data.frame(AvrBs3_deg_df$logFC)
degs_heatmap$AvrBs4 = AvrBs4_deg_df$logFC
names(degs_heatmap)[1] = "AvrBs3"

for(type in datatype){
    if(type == "pdf"){
        pdf(paste(maindir, "/images/heatmap_deg.pdf", sep=""), 5, 5, pointsize=12, pagecentre=TRUE, colormodel="rgb")
    } else {
        png(paste(maindir, "/images/heatmap_deg.png", sep=""), width=1600, height=1000, units="px", pointsize=30)
    }
    heatmap.2(as.matrix(degs_heatmap),
              main="", # heat map title
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins=c(1.8,2.2),     # widens margins around plot
              col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(9),
              dendrogram="row",
              key.xlab=expression(log[2]~fold~change),
              key.title=NA,
              adjCol=c(0.4,0.5),
              srtCol=0,
              sepwidth=c(0.01,0.01),
              sepcolor="grey",
              cexRow=1,
              cexCol=1,
              labRow = FALSE,
              ylab = "Genes",
              key.xtickfun = function() {
                  #breaks = pretty(parent.frame()$breaks)
                  #breaks = breaks[c(1, length(breaks))]
                  breaks = c(-15, -10, -5, 0, 5, 10, 15)
                  list(at = parent.frame()$scale01(breaks), labels = breaks)}
    )
    dev.off()
}

# Volcano plot
t_both_volcano = t_both_df
t_both_volcano$PValue = -log10(t_both_volcano$PValue)
deg_df_sig = deg_df[(deg_df$DEG != 0), ]
deg_df_sig$PValue = -log10(deg_df_sig$PValue)
ggplot(NULL, aes(x=logFC, y=PValue)) +
    ylab(expression(bold(-log[10]~p~value))) +
    xlab(expression(bold(log[2]~fold~change))) +
    scale_y_continuous(breaks=seq(0, 20, by=2), limits=c(0, 14.5)) +
    scale_x_continuous(breaks=seq(-10, 10, by=2), limits=c(-10, 10)) +
    plot_theme +
    geom_point(data=t_both_volcano, stat="identity", size=1, colour="grey30") +
    geom_point(data=deg_df_sig, stat="identity", size=1, colour=ggcol[2]) +
    #geom_point(data=unsig, stat="identity", size=0.6, colour="red2") +
    #geom_point(data=sig, stat="identity", size=0.6, colour="forestgreen") +
    facet_grid(. ~ treatment)
setwd(paste(maindir, "/", "images/", sep=""))
ggsave("volcano.pdf", width=20, units="cm")
ggsave("volcano.png", width=20, units="cm", dpi=450)


# Check UPA expression (https://www.ncbi.nlm.nih.gov/pubmed/19473322)
UPAs = c("LOC107875475", "LOC107847767", "LOC107856008", "LOC107867263", "LOC107864666", "LOC107858786", "LOC107871692", "LOC107864501", "LOC107875122", "LOC107869472", "LOC107866053", "LOC107857984", "UPA24", "CaBs4C.1")

data_edge_TPM$control_C = NA

UPA_rows = data_edge_TPM[which(data_edge_TPM$Gene_IDs %in% UPAs),]
UPA_rows$Gene_IDs = gsub("LOC107875475", "UPA14", UPA_rows$Gene_IDs)
UPA_rows$Gene_IDs = gsub("LOC107847767", "UPA15", UPA_rows$Gene_IDs)
UPA_rows$Gene_IDs = gsub("LOC107856008", "UPA16", UPA_rows$Gene_IDs)
UPA_rows$Gene_IDs = gsub("LOC107867263", "UPA17", UPA_rows$Gene_IDs)
UPA_rows$Gene_IDs = gsub("LOC107864666", "UPA18", UPA_rows$Gene_IDs)
UPA_rows$Gene_IDs = gsub("LOC107858786", "UPA19", UPA_rows$Gene_IDs)
UPA_rows$Gene_IDs = gsub("LOC107866053", "UPA20", UPA_rows$Gene_IDs)
UPA_rows$Gene_IDs = gsub("LOC107871692", "UPA21", UPA_rows$Gene_IDs)
UPA_rows$Gene_IDs = gsub("LOC107864501", "UPA22", UPA_rows$Gene_IDs)
UPA_rows$Gene_IDs = gsub("LOC107875122", "UPA23", UPA_rows$Gene_IDs)
UPA_rows$Gene_IDs = gsub("LOC107869472", "UPA25", UPA_rows$Gene_IDs)
UPA_rows$Gene_IDs = gsub("LOC107857984", "Bs3", UPA_rows$Gene_IDs)
UPA_rows$Gene_IDs = gsub("CaBs4C.1", "Bs4C.1", UPA_rows$Gene_IDs)
UPA_rows$Length = NULL

UPA_long = reshape(UPA_rows, varying = 2:ncol(UPA_rows), sep = "_", direction = 'long')
UPA_long_2 = melt(UPA_long)
UPA_long_2 = UPA_long_2[!(UPA_long_2$variable == "id"), ]
UPA_long_2$variable = factor(UPA_long_2$variable, levels=c("control", "AvrBs3", "AvrBs4"))


ggplot(NULL, aes(x=time, y=value, col=variable, shape=variable, order=-as.numeric(variable))) +
    geom_point(data=UPA_long_2, size=4) +
    ylab("Expression level [TPM]") +
    scale_y_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60)) +
    xlab("Replicates") +
    plot_theme +
    facet_grid(. ~ Gene_IDs) +
    scale_color_discrete(name="Conditions") +
    scale_shape_discrete(name="Conditions") +
    theme(legend.position="top")
setwd(paste(maindir, "/", "images/", sep=""))
ggsave("UPAs_Bs3.pdf", width=39, height=14, units="cm")
ggsave("UPAs_Bs3.png", width=39, units="cm", dpi=450)

gc()
reg = "up"
hitsEBE = c("LOC107838904", "LOC107840500", "LOC107840508", "LOC107843852", "LOC107847534", "LOC107857984", "LOC107862600", "LOC107862761", "LOC107866053", "LOC107869407")
hitsEBE_table = DEG_AvrBs3[DEG_AvrBs3$Gene_Symbol %in% hitsEBE, ]
pdf(paste(maindir, "/", reg, "/TSS/images/AvrBs3_", reg, "reg_predEBEs_selection.pdf", sep=""), 8, 7, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
p1 = ggplot(hitsEBE_table) +
    geom_bar(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=(2^AvrBs3_logFC)), stat="identity", colour="black", fill=ggcol[2]) +
    ylab("Fold change") +
    plot_theme +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
p2 = ggplot(hitsEBE_table) +
    geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(Control_meanTPM)), size=5, colour=ggcol[1]) +
    geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(AvrBs3_meanTPM)), size=5, colour=ggcol[2]) +
    ylab("Mean expr. level [TPM]") +
    xlab("") +
    plot_theme +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
    scale_y_continuous(limits=c(-3, 2), breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
dev.off()

hitsEBE = c("LOC107860843", "LOC107852837", "LOC107871190", "LOC107843460", "LOC107861482", "LOC107839710", "LOC107851710", "LOC107862135", "LOC107841901", "LOC107859948", "LOC107851710", "LOC107869555", "LOC107857984", "LOC107869665", "LOC107861943")
hitsEBE_table = DEG_AvrBs4[DEG_AvrBs4$Gene_Symbol %in% hitsEBE, ]
pdf(paste(maindir, "/", reg, "/TSS/images/AvrBs4_", reg, "reg_predEBEs_selection.pdf", sep=""), 10, 7, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
p1 = ggplot(hitsEBE_table) +
    geom_bar(aes(x=reorder(Gene_Symbol, -AvrBs4_logFC), y=(2^AvrBs4_logFC)), stat="identity", colour="black", fill=ggcol[3]) +
    ylab("Fold change") +
    plot_theme +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
p2 = ggplot(hitsEBE_table) +
    geom_point(aes(x=reorder(Gene_Symbol, -AvrBs4_logFC), y=-log10(Control_meanTPM)), size=5, colour=ggcol[1]) +
    geom_point(aes(x=reorder(Gene_Symbol, -AvrBs4_logFC), y=-log10(AvrBs4_meanTPM)), size=5, colour=ggcol[3]) +
    ylab("Mean expr. level [TPM]") +
    xlab("") +
    plot_theme +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
    scale_y_continuous(limits=c(-3, 2), breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
dev.off()

pdf(paste(maindir, "/", reg, "/DEG_illustrations_tables/both_", reg, "reg_FC-TPM.pdf", sep=""), 10, 7, pointsize=0.6, pagecentre=TRUE, colormodel="rgb")
#png(paste(maindir, "/DEG_illustrations_tables/both_upreg_FC-TPM.png", sep=""), width=1920, height=1080, units="px", res=230)
p1 = ggplot(test_res_both) +
    geom_bar(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=(2^AvrBs3_logFC)), stat="identity", alpha=0.8, fill=ggcol[2]) +
    geom_bar(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=(2^AvrBs4_logFC)), stat="identity", alpha=0.6, fill=ggcol[3]) +
    ylab("Fold change") +
    plot_theme +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
p2 = ggplot(test_res_both) +
    geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanControl)), size=1.5, colour=ggcol[1]) +
    geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanAvrBs3)), size=1.5, colour=ggcol[2]) +
    geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanAvrBs4)), size=1.5, colour=ggcol[3]) +
    ylab("Mean expr. level [TPM]") +
    xlab("") +
    plot_theme +
    #theme(axis.text.x=element_text(angle=90, vjust=1, hjust=1, size=8)) +
    theme(axis.text.x=element_text(angle=90, size=9, vjust=0.4)) +
    scale_y_continuous(breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
dev.off()

deg_print = test_res[test_res$AvrBs3 != 0 | test_res$AvrBs4 != 0, ]
deg_print = deg_print[-4]
deg_print = deg_print[-c(seq(6,13,1))]
deg_print[4:8] = round(deg_print[4:8], 3)
write.table(deg_print, file=paste(maindir, "/", "TALE_EBEs/deg_table_latex.txt", sep=""), col.names=FALSE, row.names=FALSE, sep=" & ", quote=FALSE)

# EBE stuff

modGFF_cmd = paste("/home/sven/Documents/Python/Projects/enrichTSS/modGFF.py -i ", genome_annot, " -o ", mod_annot, sep="")
system(modGFF_cmd)

reg_list = c("up", "down")
for(reg in reg_list){

    if(reg == "down"){

        pdf(paste(maindir, "/", reg, "/DEG_illustrations_tables/AvrBs3_", reg, "reg_FC-TPM.pdf", sep=""), 9, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
        p1 = ggplot(test_res_AvrBs3) +
            geom_bar(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=(2^AvrBs3_logFC)), stat="identity", colour="black", fill=ggcol[1]) +
            ylab("Fold change") +
            plot_theme +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
        p2 = ggplot(test_res_AvrBs3) +
            geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanControl)), size=0.5, colour=ggcol[3], position=position_jitter(width=.1, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanAvrBs3)), size=0.5, colour=ggcol[1], position=position_jitter(width=.1, height=0)) +
            ylab("Mean TPM") +
            xlab("Gene symbols") +
            plot_theme +
            theme(axis.text.x=element_text(angle=90, size=2.5, vjust=0.4)) +
            scale_y_continuous(breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
        grid.newpage()
        grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
        dev.off()

        pdf(paste(maindir, "/", reg,  "/DEG_illustrations_tables/AvrBs3_", reg, "reg_FC-TPM_meanAvrBs3.pdf", sep=""), 9, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
        p1 = ggplot(test_res_AvrBs3) +
            geom_bar(aes(x=reorder(Gene_Symbol, meanAvrBs3), y=(2^AvrBs3_logFC)), stat="identity", colour="black", fill=ggcol[1]) +
            ylab("Fold change") +
            plot_theme +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
        p2 = ggplot(test_res_AvrBs3) +
            geom_point(aes(x=reorder(Gene_Symbol, meanAvrBs3), y=-log10(meanControl)), size=0.5, colour=ggcol[3], position=position_jitter(width=.1, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, meanAvrBs3), y=-log10(meanAvrBs3)), size=0.5, colour=ggcol[1], position=position_jitter(width=.1, height=0)) +
            ylab("Mean TPM") +
            xlab("Gene symbols") +
            plot_theme +
            theme(axis.text.x=element_text(angle=90, size=2.5, vjust=0.4)) +
            scale_y_continuous(breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
        grid.newpage()
        grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
        dev.off()

        table_EBEs_AvrBs3 = NULL
        table_EBEs_AvrBs3 = AvrBs3_EBEs[,c(1, ncol(AvrBs3_EBEs), ncol(AvrBs3_EBEs)-1)]
        setDT(table_EBEs_AvrBs3)[setDT(test_res_AvrBs3), meanControl := i.meanControl, on=c("Gene_Symbol")]
        setDT(table_EBEs_AvrBs3)[setDT(test_res_AvrBs3), meanAvrBs3 := i.meanAvrBs3, on=c("Gene_Symbol")]
        table_EBEs_AvrBs3 = cbind(table_EBEs_AvrBs3, AvrBs3_EBEs[, c(2, 3, 4)])
        table_EBEs_AvrBs3 = table_EBEs_AvrBs3[with(table_EBEs_AvrBs3, order(-meanAvrBs3)), ]
        rownames(table_EBEs_AvrBs3) = NULL
        write.table(table_EBEs_AvrBs3 , file=paste(maindir, "/", reg,  "/DEG_illustrations_tables/DEG_AvrBs3_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

        table_EBEs_AvrBs3_pred = table_EBEs_AvrBs3[grep(paste(AvrBs3_EBEs_rel, collapse="|"), table_EBEs_AvrBs3$Gene_Symbol), ]
        write.table(table_EBEs_AvrBs3_pred , file=paste(maindir, "/", reg,  "/DEG_illustrations_tables/DEG_AvrBs3_EBE_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

        table_EBEs_AvrBs3_nopred = table_EBEs_AvrBs3[!grep(paste(AvrBs3_EBEs_rel, collapse="|"), table_EBEs_AvrBs3$Gene_Symbol), ]
        write.table(table_EBEs_AvrBs3_nopred , file=paste(maindir, "/", reg,  "/DEG_illustrations_tables/DEG_AvrBs3_noEBE_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)


        pdf(paste(maindir, "/", reg,  "/DEG_illustrations_tables/AvrBs4_", reg, "reg_FC-TPM.pdf", sep=""), 9, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
        p1 = ggplot(test_res_AvrBs4) +
            geom_bar(aes(x=reorder(Gene_Symbol, -AvrBs4_logFC), y=(2^AvrBs4_logFC)), stat="identity", colour="black", fill=ggcol[2]) +
            ylab("Fold change") +
            plot_theme +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
        p2 = ggplot(test_res_AvrBs4) +
            geom_point(aes(x=reorder(Gene_Symbol, -AvrBs4_logFC), y=-log10(meanControl)), size=0.5, colour=ggcol[3], position=position_jitter(width=.1, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, -AvrBs4_logFC), y=-log10(meanAvrBs4)), size=0.5, colour=ggcol[2], position=position_jitter(width=.1, height=0)) +
            ylab("Mean TPM") +
            xlab("Gene symbols") +
            plot_theme +
            theme(axis.text.x=element_text(angle=90, size=8, vjust=0.4)) +
            scale_y_continuous(breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
        grid.newpage()
        grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
        dev.off()

        pdf(paste(maindir, "/", reg, "/DEG_illustrations_tables/AvrBs4_", reg, "reg_FC-TPM_meanAvrBs4.pdf", sep=""), 9, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
        p1 = ggplot(test_res_AvrBs4) +
            geom_bar(aes(x=reorder(Gene_Symbol, meanAvrBs4), y=(2^AvrBs4_logFC)), stat="identity", colour="black", fill=ggcol[2]) +
            ylab("Fold change") +
            plot_theme +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
        p2 = ggplot(test_res_AvrBs4) +
            geom_point(aes(x=reorder(Gene_Symbol, meanAvrBs4), y=-log10(meanControl)), size=0.5, colour=ggcol[3], position=position_jitter(width=.1, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, meanAvrBs4), y=-log10(meanAvrBs4)), size=0.5, colour=ggcol[2], position=position_jitter(width=.1, height=0)) +
            ylab("Mean TPM") +
            xlab("Gene symbols") +
            plot_theme +
            theme(axis.text.x=element_text(angle=90, size=8, vjust=0.4)) +
            scale_y_continuous(breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
        grid.newpage()
        grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
        dev.off()

        AvrBs4_EBEs$Full_Name = AvrBs4_merged_full$Full_Name
        table_EBEs_AvrBs4 = AvrBs4_EBEs[,c(1, ncol(AvrBs4_EBEs), ncol(AvrBs4_EBEs)-1)]
        setDT(table_EBEs_AvrBs4)[setDT(test_res_AvrBs4), meanControl := i.meanControl, on=c("Gene_Symbol")]
        setDT(table_EBEs_AvrBs4)[setDT(test_res_AvrBs4), meanAvrBs4 := i.meanAvrBs4, on=c("Gene_Symbol")]
        table_EBEs_AvrBs4 = cbind(table_EBEs_AvrBs4, AvrBs4_EBEs[, c(2, 3, 4)])
        table_EBEs_AvrBs4 = table_EBEs_AvrBs4[with(table_EBEs_AvrBs4, order(-AvrBs4_logFC)), ]
        rownames(table_EBEs_AvrBs4) = NULL
        write.table(table_EBEs_AvrBs4 , file=paste(maindir, "/", reg, "/DEG_illustrations_tables/DEG_AvrBs4_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

        table_EBEs_AvrBs4_pred = table_EBEs_AvrBs4[grep(paste(AvrBs4_EBEs_rel, collapse="|"), table_EBEs_AvrBs4$Gene_Symbol), ]
        write.table(table_EBEs_AvrBs4_pred , file=paste(maindir, "/", reg, "/DEG_illustrations_tables/DEG_AvrBs4_EBE_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

        table_EBEs_AvrBs4_nopred = table_EBEs_AvrBs4[!grep(paste(AvrBs4_EBEs_rel, collapse="|"), table_EBEs_AvrBs4$Gene_Symbol), ]
        write.table(table_EBEs_AvrBs4_nopred , file=paste(maindir, "/", reg, "/DEG_illustrations_tables/DEG_AvrBs4_noEBE_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)


        pdf(paste(maindir, "/", reg, "/DEG_illustrations_tables/both_", reg, "reg_FC-TPM.pdf", sep=""), 9, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
        #png(paste(maindir, "/DEG_illustrations_tables/both_upreg_FC-TPM.png", sep=""), width=1920, height=1080, units="px", res=230)
        p1 = ggplot(test_res_both) +
            geom_bar(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=(2^AvrBs3_logFC)), stat="identity", alpha=0.8, fill=ggcol[1]) +
            geom_bar(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=(2^AvrBs4_logFC)), stat="identity", alpha=0.6, fill=ggcol[2]) +
            ylab("Fold change") +
            plot_theme +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
        p2 = ggplot(test_res_both) +
            geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanControl)), size=1, colour=ggcol[3], position=position_jitter(width=.2, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanAvrBs3)), size=1, colour=ggcol[1], position=position_jitter(width=.2, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanAvrBs4)), size=1, colour=ggcol[2], position=position_jitter(width=.2, height=0)) +
            ylab("Mean TPM") +
            xlab("Gene symbols") +
            plot_theme +
            theme(axis.text.x=element_text(angle=90, size=8, vjust=0.4)) +
            scale_y_continuous(breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
        grid.newpage()
        grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
        dev.off()

        pdf(paste(maindir, "/", reg, "/DEG_illustrations_tables/both_", reg, "reg_FC-TPM_Control.pdf", sep=""), 9, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
        p1 = ggplot(test_res_both) +
            geom_bar(aes(x=reorder(Gene_Symbol, meanControl), y=(2^AvrBs3_logFC)), stat="identity", alpha=0.8, fill=ggcol[1]) +
            geom_bar(aes(x=reorder(Gene_Symbol, meanControl), y=(2^AvrBs4_logFC)), stat="identity", alpha=0.6, fill=ggcol[2]) +
            ylab("Fold change") +
            plot_theme +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
        p2 = ggplot(test_res_both) +
            geom_point(aes(x=reorder(Gene_Symbol, meanControl), y=-log10(meanControl)), size=1, colour=ggcol[3], position=position_jitter(width=.2, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, meanControl), y=-log10(meanAvrBs3)), size=1, colour=ggcol[1], position=position_jitter(width=.2, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, meanControl), y=-log10(meanAvrBs4)), size=1, colour=ggcol[2], position=position_jitter(width=.2, height=0)) +
            ylab("Mean TPM") +
            xlab("Gene symbols") +
            plot_theme +
            theme(axis.text.x=element_text(angle=90, size=8, vjust=0.4)) +
            scale_y_continuous(breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
        grid.newpage()
        grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
        dev.off()

        table_EBEs_both = both_merged[, c(1, ncol(both_merged), ncol(both_merged)-2, ncol(both_merged)-1)]
        setDT(table_EBEs_both)[setDT(test_res_both), meanControl := i.meanControl, on=c("Gene_Symbol")]
        setDT(table_EBEs_both)[setDT(test_res_both), meanAvrBs3 := i.meanAvrBs3, on=c("Gene_Symbol")]
        setDT(table_EBEs_both)[setDT(test_res_both), meanAvrBs4 := i.meanAvrBs4, on=c("Gene_Symbol")]
        table_EBEs_both = cbind(table_EBEs_both, both_merged[, c(2, 3, 4, 5, 6, 7)])
        table_EBEs_both = table_EBEs_both[with(table_EBEs_both, order(-AvrBs3_logFC)), ]
        rownames(table_EBEs_both) = NULL
        write.table(table_EBEs_both , file=paste(maindir, "/", reg, "/DEG_illustrations_tables/DEG_both_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

        table_EBEs_both_pred = table_EBEs_both[grep(paste(overlap_both, collapse="|"), table_EBEs_both$Gene_Symbol), ]
        write.table(table_EBEs_both_pred , file=paste(maindir, "/", reg, "/DEG_illustrations_tables/DEG_both_EBE_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

        table_EBEs_both_nopred = table_EBEs_both[!grep(paste(overlap_both, collapse="|"), table_EBEs_both$Gene_Symbol), ]
        write.table(table_EBEs_both_nopred , file=paste(maindir, "/", reg, "/DEG_illustrations_tables/DEG_both_noEBE_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

    } else if(reg == "up"){


        EBE_home_dir = paste(maindir, "/", reg, "/", "TALE_EBEs", sep="")

        TAL_file = "AvrBs3-likes_AvrBs4-likes_TAL_file"
        file_connection = file(paste(maindir, "/", reg, "/TALE_EBEs/AvrBs3-likes_AvrBs4-likes_TAL_file", sep=""))

        AvrBs3_RVDs = "HD-NG-NS-NG-NI-NI-NI-HD-HD-NG-NS-NS-HD-HD-HD-NG-HD-NG"
        AvrBs4_RVDs = "NI-NG-NI-NI-NG-NG-NI-NS-NG-NI-NS-NG-HD-HD-NS-HD-NG-NG"
        writeLines(c(paste(">AvrBs4\t", AvrBs4_RVDs, sep=""), paste(">AvrBs3\t", AvrBs3_RVDs, sep="")), file_connection)
        close(file_connection)

        for(Avr in c("both", "AvrBs3", "AvrBs4")){

            if(Avr == "both"){
                TALVEZ_targets = 3500
                TALgetter_targets = 1500

            } else if(Avr == "AvrBs3"){
                TALVEZ_targets = 3500
                TALgetter_targets = 1500

            } else if(Avr == "AvrBs4"){
                TALVEZ_targets = 3500
                TALgetter_targets = 1500

            }

            print(Avr)

            Promotor_input = paste("ECW-123_ncbi_CIDer_v7_CDS_upstream_3000nt_", Avr, "_", reg , ".fa", sep="")
            getSeq_cmd = paste("/home/sven/Documents/Python/Projects/getSeq/getSeq.py -g ", genome_fasta, " -a ", genome_annot, " -o ", paste(EBE_home_dir, "/", Promotor_input, sep=""), " -L ", maindir, "/", reg, "/TALE_EBEs/degs_", Avr, "_", reg ,".txt -y gene -u 1250 -i 1250 -f single-line --onlyupstream", sep="")
            system(getSeq_cmd)


            TALVEZ_output_file = paste("TALVEZ_AvrBs3_AvrBs4_", Avr, "_", reg, sep="")
            TALVEZ_cmd = paste("/home/sven/Documents/Uni/Masterthesis/Analyses/TALE_tools/TALVEZ_3.2/runTALVEZ_3.2.sh -t ", TALVEZ_targets, " -l 10 -rvd ", paste(EBE_home_dir, "/", TAL_file, sep=""), " -input ", paste(EBE_home_dir, "/", Promotor_input, sep=""), " -o ", TALVEZ_output_file, sep="")
            system(TALVEZ_cmd)
            Tablemaker_output = paste(EBE_home_dir, "/results/", TALVEZ_output_file, "_tablemaker_SCORE", sep="")
            Tablemaker_cmd = paste("/home/sven/Documents/Python/Projects/EBEpred_accessories/EBEpred_tablemaker.py -p TALVEZ -c SCORE -e ", paste(EBE_home_dir, "/", TALVEZ_output_file, "_complete", sep=""), " -s pos -r 1250  -o ", Tablemaker_output, sep="")
            system(Tablemaker_cmd)


            TALgetter_output_file = paste("TALgetter_AvrBs3_AvrBs4_", Avr, "_", reg, sep="")
            TALgetter_cmd = paste("/home/sven/Documents/Python/Projects/EBEpred_accessories/runTALgetter.sh -input ", paste(EBE_home_dir, "/", Promotor_input, sep=""), "  -o ", paste(EBE_home_dir, "/", TALgetter_output_file , sep=""), " -rvd ", paste(EBE_home_dir, "/", TAL_file, sep=""), " -pthresh 0.05 -t ", TALgetter_targets, " -f true", sep="")
            system(TALgetter_cmd)
            Tablemaker_output = paste(EBE_home_dir, "/results/", TALgetter_output_file, "_tablemaker_SCORE", sep="")
            Tablemaker_cmd = paste("/home/sven/Documents/Python/Projects/EBEpred_accessories/EBEpred_tablemaker.py -p TALgetter -c Score -e ", paste(EBE_home_dir, "/", TALgetter_output_file, sep=""), " -s pos -r 1250  -o ", Tablemaker_output, sep="")
            system(Tablemaker_cmd)


            TALENT_output_file = paste("TALENT_AvrBs3_AvrBs4_", Avr, "_", reg, sep="")
            TALENT_cmd = paste("/home/sven/Documents/Python/Projects/EBEpred_accessories/runTALENT.sh -input ", paste(EBE_home_dir, "/", Promotor_input, sep=""), "  -o ", paste(EBE_home_dir, "/", TALENT_output_file , sep=""), " -rvd ", paste(EBE_home_dir, "/", TAL_file, sep=""), " -c 0 -f true -x 3.5 -n 4 -format true", sep="")
            system(TALENT_cmd)
            Tablemaker_output = paste(EBE_home_dir, "/results/", TALENT_output_file, "_tablemaker_SCORE", sep="")
            Tablemaker_cmd = paste("/home/sven/Documents/Python/Projects/EBEpred_accessories/EBEpred_tablemaker.py -p TALENT -c Score -e ", paste(EBE_home_dir, "/", TALENT_output_file, ".gff3", sep=""), " -s pos -r 1250  -o ", Tablemaker_output, sep="")
            system(Tablemaker_cmd)

            if(Avr == "both"){

                # Both analyses

                TALVEZ_both = read.table(paste(EBE_home_dir, "/results/TALVEZ_AvrBs3_AvrBs4_both_", reg,"_tablemaker_SCORE_EBEs", sep=""), header=TRUE, sep="\t")
                names(TALVEZ_both)[2:ncol(TALVEZ_both)] = paste("TALVEZ_", names(TALVEZ_both)[2:ncol(TALVEZ_both)], sep="")
                TALVEZ_both[TALVEZ_both=="n.a."] = NA

                TALgetter_both = read.table(paste(EBE_home_dir, "/results/TALgetter_AvrBs3_AvrBs4_both_", reg ,"_tablemaker_SCORE_EBEs", sep=""), header=TRUE, sep="\t")
                names(TALgetter_both)[2:ncol(TALgetter_both)] = paste("TALgetter_", names(TALgetter_both)[2:ncol(TALgetter_both)], sep="")
                TALgetter_both[TALgetter_both=="n.a."] = NA

                TALENT_both = read.table(paste(EBE_home_dir, "/results/TALENT_AvrBs3_AvrBs4_both_", reg, "_tablemaker_SCORE_EBEs", sep=""), header=TRUE, sep="\t")
                names(TALENT_both)[2:ncol(TALENT_both)] = paste("TALENT_", names(TALENT_both)[2:ncol(TALENT_both)], sep="")
                TALENT_both[TALENT_both=="n.a."] = NA

                both_merged = merge(TALVEZ_both, TALgetter_both, by="Gene_Symbol", all=TRUE)
                both_merged = merge(both_merged, TALENT_both, by="Gene_Symbol", all=TRUE)
                both_merged = merge(both_merged, AvrBs3_small_logFC, by="Gene_Symbol")
                both_merged = merge(both_merged, AvrBs4_small_logFC, by="Gene_Symbol")

                both_merged = both_merged[, c(1, grep("AvrBs3$", names(both_merged), perl=T), grep("AvrBs4$", names(both_merged), perl=T), grep("2xrep17$", names(both_merged), perl=T), grep("5\\.NS$", names(both_merged), perl=T), grep("drep9$", names(both_merged), perl=T), grep("drep15\\.16$", names(both_merged), perl=T), grep("7NS\\.drep11\\.13$", names(both_merged), perl=T), grep("drep16$", names(both_merged), perl=T), grep("5\\.NI$", names(both_merged), perl=T), grep("drep15$", names(both_merged), perl=T), grep("drep13$", names(both_merged), perl=T), grep("AvrBs3_logFC$", names(both_merged), perl=T), grep("AvrBs4_logFC$", names(both_merged), perl=T), grep("Full_Name$", names(both_merged), perl=T)) ]

                both_merged$Full_Name = mapIds(org.Cannuum.eg.db,
                                               keys=as.vector(both_merged$Gene_Symbol),
                                               column="GENENAME",
                                               keytype="SYMBOL",
                                               multiVals="first")

                write.table(both_merged , file=paste(EBE_home_dir, "/results/EBEs_both_", reg, "_all.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

                both_4 = both_merged[(is.na(both_merged$TALVEZ_AvrBs3) == 0) & (is.na(both_merged$TALgetter_AvrBs3) == 0) & (is.na(both_merged$TALVEZ_AvrBs4) == 0) & (is.na(both_merged$TALgetter_AvrBs4) == 0), ]
                write.table(both_4 , file=paste(EBE_home_dir, "/results/EBEs_both_", reg, "_4.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

                both_2 = both_merged[((is.na(both_merged$TALVEZ_AvrBs3) == 0) | (is.na(both_merged$TALgetter_AvrBs3) == 0)) & ((is.na(both_merged$TALVEZ_AvrBs4) == 0) | (is.na(both_merged$TALgetter_AvrBs4) == 0)), ]
                write.table(both_2 , file=paste(EBE_home_dir, "/results/EBEs_both_", reg, "_2.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)


            } else if(Avr == "AvrBs3"){

                # AvrBs3 Analyses

                TALVEZ_AvrBs3 = read.table(paste(EBE_home_dir, "/results/TALVEZ_AvrBs3_AvrBs4_AvrBs3_", reg, "_tablemaker_SCORE_EBEs", sep=""), header=TRUE, sep="\t")
                names(TALVEZ_AvrBs3)[2:ncol(TALVEZ_AvrBs3)] = paste("TALVEZ_", names(TALVEZ_AvrBs3)[2:ncol(TALVEZ_AvrBs3)], sep="")
                TALVEZ_AvrBs3[TALVEZ_AvrBs3=="n.a."] = NA

                TALgetter_AvrBs3 = read.table(paste(EBE_home_dir, "/results/TALgetter_AvrBs3_AvrBs4_AvrBs3_", reg, "_tablemaker_SCORE_EBEs", sep=""), header=TRUE, sep="\t")
                names(TALgetter_AvrBs3)[2:ncol(TALgetter_AvrBs3)] = paste("TALgetter_", names(TALgetter_AvrBs3)[2:ncol(TALgetter_AvrBs3)], sep="")
                TALgetter_AvrBs3[TALgetter_AvrBs3=="n.a."] = NA

                TALENT_AvrBs3 = read.table(paste(EBE_home_dir, "/results/TALENT_AvrBs3_AvrBs4_AvrBs3_", reg, "_tablemaker_SCORE_EBEs", sep=""), header=TRUE, sep="\t")
                names(TALENT_AvrBs3)[2:ncol(TALENT_AvrBs3)] = paste("TALENT_", names(TALENT_AvrBs3)[2:ncol(TALENT_AvrBs3)], sep="")
                TALENT_AvrBs3[TALENT_AvrBs3=="n.a."] = NA

                AvrBs3_merged = merge(TALVEZ_AvrBs3, TALgetter_AvrBs3, by="Gene_Symbol", all=TRUE)
                AvrBs3_merged = merge(AvrBs3_merged, TALENT_AvrBs3, by="Gene_Symbol", all=TRUE)
                AvrBs3_merged = merge(AvrBs3_merged, AvrBs3_small_logFC, by="Gene_Symbol")
                AvrBs3_merged = AvrBs3_merged[, -grep("AvrBs4" , names(AvrBs3_merged))]

                AvrBs3_merged$Full_Name = mapIds(org.Cannuum.eg.db,
                                                 keys=as.vector(AvrBs3_merged$Gene_Symbol),
                                                 column="GENENAME",
                                                 keytype="SYMBOL",
                                                 multiVals="first")
                AvrBs3_merged_full = AvrBs3_merged

                AvrBs3_merged = AvrBs3_merged[, c(1, grep("AvrBs3$", names(AvrBs3_merged), perl=T), grep("2xrep17$", names(AvrBs3_merged), perl=T), grep("5\\.NS$", names(AvrBs3_merged), perl=T), grep("drep9$", names(AvrBs3_merged), perl=T), grep("drep15\\.16$", names(AvrBs3_merged), perl=T), grep("7NS\\.drep11\\.13$", names(AvrBs3_merged), perl=T), grep("drep16$", names(AvrBs3_merged), perl=T), grep("AvrBs3_logFC$", names(AvrBs3_merged), perl=T), grep("Full_Name$", names(AvrBs3_merged), perl=T)) ]
                AvrBs3_EBEs = AvrBs3_merged
                write.table(AvrBs3_merged , file=paste(EBE_home_dir, "/results/EBEs_AvrBs3_", reg, "_all.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

                AvrBs3_3 = AvrBs3_merged[(is.na(AvrBs3_merged$TALVEZ_AvrBs3) == 0) & (is.na(AvrBs3_merged$TALgetter_AvrBs3) == 0) & (is.na(AvrBs3_merged$TALENT_AvrBs3) == 0), ]
                write.table(AvrBs3_3 , file=paste(EBE_home_dir, "/results/EBEs_AvrBs3_", reg, "_3.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

                AvrBs3_TALVEZ_TALgetter = AvrBs3_merged[(is.na(AvrBs3_merged$TALVEZ_AvrBs3) == 0) & (is.na(AvrBs3_merged$TALgetter_AvrBs3) == 0) & (is.na(AvrBs3_merged$TALENT_AvrBs3) == 1), ]
                AvrBs3_TALVEZ_TALENT = AvrBs3_merged[(is.na(AvrBs3_merged$TALVEZ_AvrBs3) == 0) & (is.na(AvrBs3_merged$TALgetter_AvrBs3) == 1) & (is.na(AvrBs3_merged$TALENT_AvrBs3) == 0), ]
                AvrBs3_TALgetter_TALENT = AvrBs3_merged[(is.na(AvrBs3_merged$TALVEZ_AvrBs3) == 1) & (is.na(AvrBs3_merged$TALgetter_AvrBs3) == 0) & (is.na(AvrBs3_merged$TALENT_AvrBs3) == 0), ]
                AvrBs3_2 = rbind(AvrBs3_TALVEZ_TALgetter, AvrBs3_TALVEZ_TALENT)
                AvrBs3_2 = rbind(AvrBs3_2, AvrBs3_TALgetter_TALENT)
                write.table(AvrBs3_2 , file=paste(EBE_home_dir, "/results/EBEs_AvrBs3_", reg, "_2.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

                AvrBs3_none = AvrBs3_merged[(is.na(AvrBs3_merged$TALVEZ_AvrBs3) == 1) & (is.na(AvrBs3_merged$TALgetter_AvrBs3) == 1) & (is.na(AvrBs3_merged$TALENT_AvrBs3) == 1), ]
                write.table(AvrBs3_none , file=paste(EBE_home_dir, "/results/EBEs_AvrBs3_", reg, "_none.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

                AvrBs3_TALVEZ = AvrBs3_merged[(is.na(AvrBs3_merged$TALVEZ_AvrBs3) == 0) & (is.na(AvrBs3_merged$TALgetter_AvrBs3) == 1) & (is.na(AvrBs3_merged$TALENT_AvrBs3) == 1), ]
                AvrBs3_TALgetter = AvrBs3_merged[(is.na(AvrBs3_merged$TALVEZ_AvrBs3) == 1) & (is.na(AvrBs3_merged$TALgetter_AvrBs3) == 0) & (is.na(AvrBs3_merged$TALENT_AvrBs3) == 1), ]
                AvrBs3_TALENT = AvrBs3_merged[(is.na(AvrBs3_merged$TALVEZ_AvrBs3) == 1) & (is.na(AvrBs3_merged$TALgetter_AvrBs3) == 1) & (is.na(AvrBs3_merged$TALENT_AvrBs3) == 0), ]
                AvrBs3_1 = rbind(AvrBs3_none, AvrBs3_TALVEZ)
                AvrBs3_1 = rbind(AvrBs3_1, AvrBs3_TALgetter)
                AvrBs3_1 = rbind(AvrBs3_1, AvrBs3_TALENT)
                write.table(AvrBs3_1 , file=paste(EBE_home_dir, "/results/EBEs_AvrBs3_", reg, "_01.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)


            } else if(Avr == "AvrBs4"){

                # AvrBs4 Analyses

                TALVEZ_AvrBs4 = read.table(paste(EBE_home_dir, "/results/TALVEZ_AvrBs3_AvrBs4_AvrBs4_", reg, "_tablemaker_SCORE_EBEs", sep=""), header=TRUE, sep="\t")
                names(TALVEZ_AvrBs4)[2:ncol(TALVEZ_AvrBs4)] = paste("TALVEZ_", names(TALVEZ_AvrBs4)[2:ncol(TALVEZ_AvrBs4)], sep="")
                TALVEZ_AvrBs4[TALVEZ_AvrBs4=="n.a."] = NA

                TALgetter_AvrBs4 = read.table(paste(EBE_home_dir, "/results/TALgetter_AvrBs3_AvrBs4_AvrBs4_", reg, "_tablemaker_SCORE_EBEs", sep=""), header=TRUE, sep="\t")
                names(TALgetter_AvrBs4)[2:ncol(TALgetter_AvrBs4)] = paste("TALgetter_", names(TALgetter_AvrBs4)[2:ncol(TALgetter_AvrBs4)], sep="")
                TALgetter_AvrBs4[TALgetter_AvrBs4=="n.a."] = NA

                TALENT_AvrBs4 = read.table(paste(EBE_home_dir, "/results/TALENT_AvrBs3_AvrBs4_AvrBs4_", reg, "_tablemaker_SCORE_EBEs", sep=""), header=TRUE, sep="\t")
                names(TALENT_AvrBs4)[2:ncol(TALENT_AvrBs4)] = paste("TALENT_", names(TALENT_AvrBs4)[2:ncol(TALENT_AvrBs4)], sep="")
                TALENT_AvrBs4[TALENT_AvrBs4=="n.a."] = NA

                AvrBs4_merged = merge(TALVEZ_AvrBs4, TALgetter_AvrBs4, by="Gene_Symbol", all=TRUE)
                AvrBs4_merged = merge(AvrBs4_merged, TALENT_AvrBs4, by="Gene_Symbol", all=TRUE)
                AvrBs4_merged = merge(AvrBs4_merged, AvrBs4_small_logFC, by="Gene_Symbol")
                AvrBs4_merged = AvrBs4_merged[, -grep("AvrBs3" , names(AvrBs4_merged))]

                AvrBs4_merged$Full_Name = mapIds(org.Cannuum.eg.db,
                                                 keys=as.vector(AvrBs4_merged$Gene_Symbol),
                                                 column="GENENAME",
                                                 keytype="SYMBOL",
                                                 multiVals="first")
                AvrBs4_merged_full = AvrBs4_merged

                AvrBs4_merged = AvrBs4_merged[, c(1, grep("AvrBs4$", names(AvrBs4_merged), perl=T), grep("5\\.NI$", names(AvrBs4_merged), perl=T), grep("drep15$", names(AvrBs4_merged), perl=T), grep("drep13$", names(AvrBs4_merged), perl=T), grep("AvrBs4_logFC$", names(AvrBs4_merged), perl=T), grep("Full_Name$", names(AvrBs4_merged), perl=T)) ]
                AvrBs4_EBEs = AvrBs4_merged
                write.table(AvrBs4_merged , file=paste(EBE_home_dir, "/results/EBEs_AvrBs4_", reg, "_all.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

                AvrBs4_3 = AvrBs4_merged[(is.na(AvrBs4_merged$TALVEZ_AvrBs4) == 0) & (is.na(AvrBs4_merged$TALgetter_AvrBs4) == 0) & (is.na(AvrBs4_merged$TALENT_AvrBs4) == 0), ]
                write.table(AvrBs4_3 , file=paste(EBE_home_dir, "/results/EBEs_AvrBs4_", reg, "_3.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

                AvrBs4_TALVEZ_TALgetter = AvrBs4_merged[(is.na(AvrBs4_merged$TALVEZ_AvrBs4) == 0) & (is.na(AvrBs4_merged$TALgetter_AvrBs4) == 0) & (is.na(AvrBs4_merged$TALENT_AvrBs4) == 1), ]
                AvrBs4_TALVEZ_TALENT = AvrBs4_merged[(is.na(AvrBs4_merged$TALVEZ_AvrBs4) == 0) & (is.na(AvrBs4_merged$TALgetter_AvrBs4) == 1) & (is.na(AvrBs4_merged$TALENT_AvrBs4) == 0), ]
                AvrBs4_TALgetter_TALENT = AvrBs4_merged[(is.na(AvrBs4_merged$TALVEZ_AvrBs4) == 1) & (is.na(AvrBs4_merged$TALgetter_AvrBs4) == 0) & (is.na(AvrBs4_merged$TALENT_AvrBs4) == 0), ]
                AvrBs4_2 = rbind(AvrBs4_TALVEZ_TALgetter, AvrBs4_TALVEZ_TALENT)
                AvrBs4_2 = rbind(AvrBs4_2, AvrBs4_TALgetter_TALENT)
                write.table(AvrBs4_2 , file=paste(EBE_home_dir, "/results/EBEs_AvrBs4_", reg, "_2.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

                AvrBs4_none = AvrBs4_merged[(is.na(AvrBs4_merged$TALVEZ_AvrBs4) == 1) & (is.na(AvrBs4_merged$TALgetter_AvrBs4) == 1), ]
                write.table(AvrBs4_none , file=paste(EBE_home_dir, "/results/EBEs_AvrBs4_", reg, "_none.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

                AvrBs4_TALVEZ = AvrBs4_merged[(is.na(AvrBs4_merged$TALVEZ_AvrBs4) == 0) & (is.na(AvrBs4_merged$TALgetter_AvrBs4) == 1) & (is.na(AvrBs4_merged$TALENT_AvrBs4) == 1), ]
                AvrBs4_TALgetter = AvrBs4_merged[(is.na(AvrBs4_merged$TALVEZ_AvrBs4) == 1) & (is.na(AvrBs4_merged$TALgetter_AvrBs4) == 0) & (is.na(AvrBs4_merged$TALENT_AvrBs4) == 1), ]
                AvrBs4_TALENT = AvrBs4_merged[(is.na(AvrBs4_merged$TALVEZ_AvrBs4) == 1) & (is.na(AvrBs4_merged$TALgetter_AvrBs4) == 1) & (is.na(AvrBs4_merged$TALENT_AvrBs4) == 0), ]
                AvrBs4_1 = rbind(AvrBs4_none, AvrBs4_TALVEZ)
                AvrBs4_1 = rbind(AvrBs4_1, AvrBs4_TALgetter)
                AvrBs4_1 = rbind(AvrBs4_1, AvrBs4_TALENT)
                write.table(AvrBs4_1 , file=paste(EBE_home_dir, "/results/EBEs_AvrBs4_", reg, "_01.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
            }

            pwr_mat = read.table(file="/home/sven/Documents/Uni/Masterthesis/Analyses/binding_affinity_rvds.tsv", sep="\t", header=TRUE, na.strings="n.a.")
            UPA20_EBE = "ATATAAACCTNNCCCTCT"
            CaBs4C1_EBE = "ATAATTANTANTCCNCTT"

            if(Avr != "both"){
                TALVEZ_EBEs = readDNAStringSet(paste(EBE_home_dir, "/results/EBE_sequences/", TALVEZ_output_file, "_tablemaker_SCORE_EBEs_", Avr, ".fa", sep=""))
                TALgetter_EBEs = readDNAStringSet(paste(EBE_home_dir, "/results/EBE_sequences/", TALgetter_output_file, "_tablemaker_SCORE_EBEs_", Avr, ".fa", sep=""))
                TALENT_EBEs = readDNAStringSet(paste(EBE_home_dir, "/results/EBE_sequences/", TALENT_output_file, "_tablemaker_SCORE_EBEs_", Avr, ".fa", sep=""))

                #TALE = readDNAStringSet(paste("/home/sven/Documents/Uni/Masterthesis/Analyses/", Avr, ".fa", sep=""))

                Avr_EBEs = append(TALVEZ_EBEs, TALgetter_EBEs)
                Avr_EBEs = append(Avr_EBEs, TALENT_EBEs)

                names(Avr_EBEs) = sub("((AvrBs3)|(AvrBs4))_", "", names(Avr_EBEs))
                names(Avr_EBEs) = sub("((TALENT)|(TALVEZ)|(TALgetter))_", "", names(Avr_EBEs))
                names(Avr_EBEs) = sub("\\.cds.+?\\.1", "", names(Avr_EBEs))
                names(Avr_EBEs) = sub("\\.gene.+?\\.1", "", names(Avr_EBEs))
                Avr_EBEs = Avr_EBEs[!duplicated(names(Avr_EBEs))]
                Avr_EBEs = DNAStringSet(Avr_EBEs, start=2)

                if(Avr == "AvrBs3"){
                    rel_genes = rbind(AvrBs3_3, AvrBs3_2, AvrBs3_1)
                    Avr_EBEs = Avr_EBEs[grep(paste(as.vector(rel_genes$Gene_Symbol), collapse="|"), names(Avr_EBEs), perl=TRUE)]
                } else if(Avr == "AvrBs4"){
                    rel_genes = rbind(AvrBs4_3, AvrBs4_2, AvrBs4_1)
                    Avr_EBEs = Avr_EBEs[grep(paste(as.vector(rel_genes$Gene_Symbol), collapse="|"), names(Avr_EBEs), perl=TRUE)]
                }

                ii = 1
                keep = NULL
                for(pos in 1:length(Avr_EBEs)){
                    sequence = Avr_EBEs[pos]
                    if(Avr == "AvrBs3"){
                        result_evalEBE = eval_EBE(RVD_set=AvrBs3_RVDs, EBE=as.character(sequence[[1]]), pwr_mat, relative_EBE=UPA20_EBE, first_x_perfect=4)
                        if(result_evalEBE[2] >= 2 & result_evalEBE[3] <= 4 & result_evalEBE[4] == 0){
                            keep[ii] = pos
                            ii = ii + 1
                        }

                        # if(grepl("(A|C)T(A|C|G)T", as.character(sequence[[1]][1:5]), perl=TRUE)){
                        #     keep[ii] = pos
                        #     ii = ii + 1
                        # }

                    } else if(Avr == "AvrBs4"){
                        result_evalEBE = eval_EBE(RVD_set=AvrBs4_RVDs, EBE=as.character(sequence[[1]]), pwr_mat, relative_EBE=CaBs4C1_EBE, first_x_perfect=2)
                        if(result_evalEBE[2] >= 2 & result_evalEBE[3] <= 4 & result_evalEBE[4] == 0){
                            keep[ii] = pos
                            ii = ii + 1
                        }


                        # if(grepl("AT", as.character(sequence[[1]][1:2]), perl=TRUE)){
                        #     keep[ii] = pos
                        #     ii = ii + 1
                        # }
                    }
                }

                Avr_EBEs = Avr_EBEs[keep]

                #Avr_EBEs = append(TALE, Avr_EBEs)

                # alignment = msa(Avr_EBEs, method="ClustalW", gapOpening=100, maxiters=25, order="input")
                # setwd(paste(EBE_home_dir, "/results/EBE_alignments", sep=""))
                # msaPrettyPrint(alignment, file=paste(Avr,"_", reg, ".tex", sep="") , output="tex", showLogo="top", paperHeight=9, verbose=FALSE, askForOverwrite=FALSE, consensusColor="ColdHot", shadingMode="similar")
                #
                # seqinr_alignment = msaConvert(alignment, type="seqinr::alignment")
                # d = dist.alignment(seqinr_alignment, "identity")
                # EBE_tree = nj(d)
                #
                # pdf(paste(EBE_home_dir, "/results/EBE_alignments/", Avr, "_", reg, "_tree.pdf", sep=""), 5, 8, pointsize=9, pagecentre=TRUE, colormodel="rgb")
                # plot(EBE_tree)
                # dev.off()

                cleannames_Avr_EBEs = sub("\\..+", "", names(Avr_EBEs))
                cleannames_Avr_EBEs = cleannames_Avr_EBEs[2:length(cleannames_Avr_EBEs)]
                write.table(cleannames_Avr_EBEs, file=paste(EBE_home_dir, "/results/EBE_alignments/", Avr, "_", reg, "_EBEs_relevant.txt", sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

                if(Avr == "AvrBs3"){
                    AvrBs3_EBEs_rel = unique(cleannames_Avr_EBEs)
                    AvrBs3_EBEs_orig = Avr_EBEs
                } else if(Avr == "AvrBs4"){
                    AvrBs4_EBEs_rel = unique(cleannames_Avr_EBEs)
                    AvrBs4_EBEs_orig = Avr_EBEs
                }


            } else if(Avr == "both"){
                GoI_both = c(as.vector(both_4$Gene_Symbol), as.vector(both_2$Gene_Symbol))

                TALVEZ_EBEs_AvrBs3 = readDNAStringSet(paste(EBE_home_dir, "/results/EBE_sequences/", TALVEZ_output_file, "_tablemaker_SCORE_EBEs_AvrBs3.fa", sep=""))
                TALgetter_EBEss_AvrBs3 = readDNAStringSet(paste(EBE_home_dir, "/results/EBE_sequences/", TALgetter_output_file, "_tablemaker_SCORE_EBEs_AvrBs3.fa", sep=""))
                TALENT_EBEss_AvrBs3 = readDNAStringSet(paste(EBE_home_dir, "/results/EBE_sequences/", TALENT_output_file, "_tablemaker_SCORE_EBEs_AvrBs3.fa", sep=""))

                #TALE_AvrBs3 = readDNAStringSet(paste("/home/sven/Documents/Uni/Masterthesis/Analyses/AvrBs3.fa", sep=""))

                AvrBs3_both_EBEs = append(TALVEZ_EBEs_AvrBs3, TALgetter_EBEss_AvrBs3)
                AvrBs3_both_EBEs = append(AvrBs3_both_EBEs, TALENT_EBEss_AvrBs3)

                AvrBs3_both_EBEs = AvrBs3_both_EBEs[unique(grep(paste(GoI_both, collapse="|"), names(AvrBs3_both_EBEs)))]
                names(AvrBs3_both_EBEs) = sub("AvrBs3_", "", names(AvrBs3_both_EBEs))
                names(AvrBs3_both_EBEs) = sub("((TALENT)|(TALVEZ)|(TALgetter))_", "", names(AvrBs3_both_EBEs))
                names(AvrBs3_both_EBEs) = sub("\\.cds.+?\\.1", "", names(AvrBs3_both_EBEs))
                names(AvrBs3_both_EBEs) = sub("\\.gene.+?\\.1", "", names(AvrBs3_both_EBEs))
                AvrBs3_both_EBEs = AvrBs3_both_EBEs[!duplicated(names(AvrBs3_both_EBEs))]
                AvrBs3_both_EBEs = DNAStringSet(AvrBs3_both_EBEs, start=2)

                ii = 1
                keep = NULL
                for(pos in 1:length(AvrBs3_both_EBEs)){
                    sequence = AvrBs3_both_EBEs[pos]
                    result_evalEBE = eval_EBE(RVD_set=AvrBs3_RVDs, EBE=as.character(sequence[[1]]), pwr_mat, relative_EBE=UPA20_EBE, first_x_perfect=4)
                    if(result_evalEBE[2] >= 2 & result_evalEBE[3] <= 4 & result_evalEBE[4] == 0){
                        keep[ii] = pos
                        ii = ii + 1
                    }

                    # if(grepl("(A|C)T(A|C|G)T", as.character(sequence[[1]][1:5]), perl=TRUE)){
                    #     keep[ii] = pos
                    #     ii = ii + 1
                    # }
                }

                AvrBs3_both_EBEs = AvrBs3_both_EBEs[keep]
                #AvrBs3_both_EBEs_Avr = append(TALE_AvrBs3, AvrBs3_both_EBEs)

                #alignment = msa(AvrBs3_both_EBEs_Avr, method="ClustalW", gapOpening=100, maxiters=25, order="input")
                #setwd(paste(EBE_home_dir, "/results/", "EBE_alignments", sep=""))
                #msaPrettyPrint(alignment, file=paste("both_", reg, "_AvrBs3.tex", sep="") , output="tex", showLogo="top", paperHeight=6, verbose=FALSE, askForOverwrite=FALSE, consensusColor="ColdHot", shadingMode="similar")

                # seqinr_alignment = msaConvert(alignment, type="seqinr::alignment")
                # d = dist.alignment(seqinr_alignment, "identity")
                # EBE_tree = nj(d)
                #
                # pdf(paste(EBE_home_dir, "/results/EBE_alignments/both_", reg, "_AvrBs3_tree.pdf", sep=""), 5, 5, pointsize=9, pagecentre=TRUE, colormodel="rgb")
                # plot(EBE_tree)
                # dev.off()


                TALVEZ_EBEs_AvrBs4 = readDNAStringSet(paste(EBE_home_dir, "/results/EBE_sequences/", TALVEZ_output_file, "_tablemaker_SCORE_EBEs_AvrBs4.fa", sep=""))
                TALgetter_EBEss_AvrBs4 = readDNAStringSet(paste(EBE_home_dir, "/results/EBE_sequences/", TALgetter_output_file, "_tablemaker_SCORE_EBEs_AvrBs4.fa", sep=""))
                TALENT_EBEss_AvrBs4 = readDNAStringSet(paste(EBE_home_dir, "/results/EBE_sequences/", TALENT_output_file, "_tablemaker_SCORE_EBEs_AvrBs4.fa", sep=""))

                #TALE_AvrBs4 = readDNAStringSet(paste("/home/sven/Documents/Uni/Masterthesis/Analyses/AvrBs4.fa", sep=""))

                AvrBs4_both_EBEs = append(TALVEZ_EBEs_AvrBs4, TALgetter_EBEss_AvrBs4)
                AvrBs4_both_EBEs = append(AvrBs4_both_EBEs, TALENT_EBEss_AvrBs4)

                AvrBs4_both_EBEs = AvrBs4_both_EBEs[unique(grep(paste(GoI_both, collapse="|"), names(AvrBs4_both_EBEs)))]
                names(AvrBs4_both_EBEs) = sub("AvrBs4_", "", names(AvrBs4_both_EBEs))
                names(AvrBs4_both_EBEs) = sub("((TALENT)|(TALVEZ)|(TALgetter))_", "", names(AvrBs4_both_EBEs))
                names(AvrBs4_both_EBEs) = sub("\\.cds.+?\\.1", "", names(AvrBs4_both_EBEs))
                names(AvrBs4_both_EBEs) = sub("\\.gene.+?\\.1", "", names(AvrBs4_both_EBEs))
                AvrBs4_both_EBEs = AvrBs4_both_EBEs[!duplicated(names(AvrBs4_both_EBEs))]
                AvrBs4_both_EBEs = DNAStringSet(AvrBs4_both_EBEs, start=2)

                ii = 1
                keep = NULL
                for(pos in 1:length(AvrBs4_both_EBEs)){
                    sequence = AvrBs4_both_EBEs[pos]
                    result_evalEBE = eval_EBE(RVD_set=AvrBs4_RVDs, EBE=as.character(sequence[[1]]), pwr_mat, relative_EBE=CaBs4C1_EBE, first_x_perfect=2)
                    if(result_evalEBE[2] >= 2 & result_evalEBE[3] <= 4 & result_evalEBE[4] == 0){
                        keep[ii] = pos
                        ii = ii + 1
                    }

                    # if(grepl("AT", as.character(sequence[[1]][1:2]), perl=TRUE)){
                    #     keep[ii] = pos
                    #     ii = ii + 1
                    # }
                }

                AvrBs4_both_EBEs = AvrBs4_both_EBEs[keep]
                #AvrBs4_both_EBEs_Avr = append(TALE_AvrBs4, AvrBs4_both_EBEs)

                #alignment = msa(AvrBs4_both_EBEs_Avr, method="ClustalW", gapOpening=100, maxiters=25, order="input")
                #setwd(paste(EBE_home_dir, "/results/", "EBE_alignments", sep=""))
                #msaPrettyPrint(alignment, file=paste("both_AvrBs4.tex", sep="") , output="tex", showLogo="top", paperHeight=4, verbose=FALSE, askForOverwrite=FALSE, consensusColor="ColdHot", shadingMode="similar")
                #
                #seqinr_alignment = msaConvert(alignment, type="seqinr::alignment")
                #d = dist.alignment(seqinr_alignment, "identity")
                #EBE_tree = nj(d)
                #
                #pdf(paste(EBE_home_dir, "/results/EBE_alignments/both_", reg, "_AvrBs4_tree.pdf", sep=""), 5, 5, pointsize=9, pagecentre=TRUE, colormodel="rgb")
                #plot(EBE_tree)
                #dev.off()

                both_EBEs_orig_AvrBs3 = AvrBs3_both_EBEs
                both_EBEs_orig_AvrBs4 = AvrBs4_both_EBEs

                overlap_ready_AvrBs3 = names(AvrBs3_both_EBEs)
                overlap_ready_AvrBs3 = sub("\\..+", "", overlap_ready_AvrBs3)
                overlap_ready_AvrBs4 = names(AvrBs4_both_EBEs)
                overlap_ready_AvrBs4 = sub("\\..+", "", overlap_ready_AvrBs4)

                overlap_both = c(unique(overlap_ready_AvrBs3), unique(overlap_ready_AvrBs4))
                overlap_both = sub("\\..+", "", overlap_both)
                overlap_both = unique(overlap_both[duplicated(overlap_both)])

                write.table(overlap_both , file=paste(EBE_home_dir, "/results/EBE_alignments/both_", reg, "_EBEs_overlap.txt", sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

            }
        }


        pdf(paste(maindir, "/", reg, "/DEG_illustrations_tables/AvrBs3_", reg, "reg_FC-TPM.pdf", sep=""), 9, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
        p1 = ggplot(test_res_AvrBs3) +
            geom_bar(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=(2^AvrBs3_logFC)), stat="identity", colour="black", fill=ggcol[1]) +
            ylab("Fold change") +
            plot_theme +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
        p2 = ggplot(test_res_AvrBs3) +
            geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanControl)), size=0.5, colour=ggcol[3], position=position_jitter(width=.1, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanAvrBs3)), size=0.5, colour=ggcol[1], position=position_jitter(width=.1, height=0)) +
            ylab("Mean TPM") +
            xlab("Gene symbols") +
            plot_theme +
            theme(axis.text.x=element_text(angle=90, size=2.5, vjust=0.4)) +
            scale_y_continuous(breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
        grid.newpage()
        grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
        dev.off()

        pdf(paste(maindir, "/", reg,  "/DEG_illustrations_tables/AvrBs3_", reg, "reg_FC-TPM_meanAvrBs3.pdf", sep=""), 9, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
        p1 = ggplot(test_res_AvrBs3) +
            geom_bar(aes(x=reorder(Gene_Symbol, meanAvrBs3), y=(2^AvrBs3_logFC)), stat="identity", colour="black", fill=ggcol[1]) +
            ylab("Fold change") +
            plot_theme +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
        p2 = ggplot(test_res_AvrBs3) +
            geom_point(aes(x=reorder(Gene_Symbol, meanAvrBs3), y=-log10(meanControl)), size=0.5, colour=ggcol[3], position=position_jitter(width=.1, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, meanAvrBs3), y=-log10(meanAvrBs3)), size=0.5, colour=ggcol[1], position=position_jitter(width=.1, height=0)) +
            ylab("Mean TPM") +
            xlab("Gene symbols") +
            plot_theme +
            theme(axis.text.x=element_text(angle=90, size=2.5, vjust=0.4)) +
            scale_y_continuous(breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
        grid.newpage()
        grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
        dev.off()

        table_EBEs_AvrBs3 = NULL
        table_EBEs_AvrBs3 = AvrBs3_EBEs[,c(1, ncol(AvrBs3_EBEs), ncol(AvrBs3_EBEs)-1)]
        setDT(table_EBEs_AvrBs3)[setDT(test_res_AvrBs3), meanControl := i.meanControl, on=c("Gene_Symbol")]
        setDT(table_EBEs_AvrBs3)[setDT(test_res_AvrBs3), meanAvrBs3 := i.meanAvrBs3, on=c("Gene_Symbol")]
        table_EBEs_AvrBs3 = cbind(table_EBEs_AvrBs3, AvrBs3_EBEs[, c(2, 3, 4)])
        table_EBEs_AvrBs3 = table_EBEs_AvrBs3[with(table_EBEs_AvrBs3, order(-meanAvrBs3)), ]
        rownames(table_EBEs_AvrBs3) = NULL
        write.table(table_EBEs_AvrBs3 , file=paste(maindir, "/", reg,  "/DEG_illustrations_tables/DEG_AvrBs3_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

        table_EBEs_AvrBs3_pred = table_EBEs_AvrBs3[grep(paste(AvrBs3_EBEs_rel, collapse="|"), table_EBEs_AvrBs3$Gene_Symbol), ]
        write.table(table_EBEs_AvrBs3_pred , file=paste(maindir, "/", reg,  "/DEG_illustrations_tables/DEG_AvrBs3_EBE_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

        table_EBEs_AvrBs3_nopred = table_EBEs_AvrBs3[!grep(paste(AvrBs3_EBEs_rel, collapse="|"), table_EBEs_AvrBs3$Gene_Symbol), ]
        write.table(table_EBEs_AvrBs3_nopred , file=paste(maindir, "/", reg,  "/DEG_illustrations_tables/DEG_AvrBs3_noEBE_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)



        pdf(paste(maindir, "/", reg,  "/DEG_illustrations_tables/AvrBs4_", reg, "reg_FC-TPM.pdf", sep=""), 9, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
        p1 = ggplot(test_res_AvrBs4) +
            geom_bar(aes(x=reorder(Gene_Symbol, -AvrBs4_logFC), y=(2^AvrBs4_logFC)), stat="identity", colour="black", fill=ggcol[2]) +
            ylab("Fold change") +
            plot_theme +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
        p2 = ggplot(test_res_AvrBs4) +
            geom_point(aes(x=reorder(Gene_Symbol, -AvrBs4_logFC), y=-log10(meanControl)), size=0.5, colour=ggcol[3], position=position_jitter(width=.1, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, -AvrBs4_logFC), y=-log10(meanAvrBs4)), size=0.5, colour=ggcol[2], position=position_jitter(width=.1, height=0)) +
            ylab("Mean TPM") +
            xlab("Gene symbols") +
            plot_theme +
            theme(axis.text.x=element_text(angle=90, size=8, vjust=0.4)) +
            scale_y_continuous(breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
        grid.newpage()
        grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
        dev.off()

        pdf(paste(maindir, "/", reg, "/DEG_illustrations_tables/AvrBs4_", reg, "reg_FC-TPM_meanAvrBs4.pdf", sep=""), 9, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
        p1 = ggplot(test_res_AvrBs4) +
            geom_bar(aes(x=reorder(Gene_Symbol, meanAvrBs4), y=(2^AvrBs4_logFC)), stat="identity", colour="black", fill=ggcol[2]) +
            ylab("Fold change") +
            plot_theme +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
        p2 = ggplot(test_res_AvrBs4) +
            geom_point(aes(x=reorder(Gene_Symbol, meanAvrBs4), y=-log10(meanControl)), size=0.5, colour=ggcol[3], position=position_jitter(width=.1, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, meanAvrBs4), y=-log10(meanAvrBs4)), size=0.5, colour=ggcol[2], position=position_jitter(width=.1, height=0)) +
            ylab("Mean TPM") +
            xlab("Gene symbols") +
            plot_theme +
            theme(axis.text.x=element_text(angle=90, size=8, vjust=0.4)) +
            scale_y_continuous(breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
        grid.newpage()
        grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
        dev.off()

        AvrBs4_EBEs$Full_Name = AvrBs4_merged_full$Full_Name
        table_EBEs_AvrBs4 = AvrBs4_EBEs[,c(1, ncol(AvrBs4_EBEs), ncol(AvrBs4_EBEs)-1)]
        setDT(table_EBEs_AvrBs4)[setDT(test_res_AvrBs4), meanControl := i.meanControl, on=c("Gene_Symbol")]
        setDT(table_EBEs_AvrBs4)[setDT(test_res_AvrBs4), meanAvrBs4 := i.meanAvrBs4, on=c("Gene_Symbol")]
        table_EBEs_AvrBs4 = cbind(table_EBEs_AvrBs4, AvrBs4_EBEs[, c(2, 3, 4)])
        table_EBEs_AvrBs4 = table_EBEs_AvrBs4[with(table_EBEs_AvrBs4, order(-AvrBs4_logFC)), ]
        rownames(table_EBEs_AvrBs4) = NULL
        write.table(table_EBEs_AvrBs4 , file=paste(maindir, "/", reg, "/DEG_illustrations_tables/DEG_AvrBs4_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

        table_EBEs_AvrBs4_pred = table_EBEs_AvrBs4[grep(paste(AvrBs4_EBEs_rel, collapse="|"), table_EBEs_AvrBs4$Gene_Symbol), ]
        write.table(table_EBEs_AvrBs4_pred , file=paste(maindir, "/", reg, "/DEG_illustrations_tables/DEG_AvrBs4_EBE_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

        table_EBEs_AvrBs4_nopred = table_EBEs_AvrBs4[!grep(paste(AvrBs4_EBEs_rel, collapse="|"), table_EBEs_AvrBs4$Gene_Symbol), ]
        write.table(table_EBEs_AvrBs4_nopred , file=paste(maindir, "/", reg, "/DEG_illustrations_tables/DEG_AvrBs4_noEBE_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)



        pdf(paste(maindir, "/", reg, "/DEG_illustrations_tables/both_", reg, "reg_FC-TPM.pdf", sep=""), 9, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
        #png(paste(maindir, "/DEG_illustrations_tables/both_upreg_FC-TPM.png", sep=""), width=1920, height=1080, units="px", res=230)
        p1 = ggplot(test_res_both) +
            geom_bar(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=(2^AvrBs3_logFC)), stat="identity", alpha=0.8, fill=ggcol[1]) +
            geom_bar(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=(2^AvrBs4_logFC)), stat="identity", alpha=0.6, fill=ggcol[2]) +
            ylab("Fold change") +
            plot_theme +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
        p2 = ggplot(test_res_both) +
            geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanControl)), size=1, colour=ggcol[3], position=position_jitter(width=.2, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanAvrBs3)), size=1, colour=ggcol[1], position=position_jitter(width=.2, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanAvrBs4)), size=1, colour=ggcol[2], position=position_jitter(width=.2, height=0)) +
            ylab("Mean TPM") +
            xlab("Gene symbols") +
            plot_theme +
            theme(axis.text.x=element_text(angle=90, size=8, vjust=0.4)) +
            scale_y_continuous(breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
        grid.newpage()
        grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
        dev.off()

        pdf(paste(maindir, "/", reg, "/DEG_illustrations_tables/both_", reg, "reg_FC-TPM_Control.pdf", sep=""), 9, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
        p1 = ggplot(test_res_both) +
            geom_bar(aes(x=reorder(Gene_Symbol, meanControl), y=(2^AvrBs3_logFC)), stat="identity", alpha=0.8, fill=ggcol[1]) +
            geom_bar(aes(x=reorder(Gene_Symbol, meanControl), y=(2^AvrBs4_logFC)), stat="identity", alpha=0.6, fill=ggcol[2]) +
            ylab("Fold change") +
            plot_theme +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
        p2 = ggplot(test_res_both) +
            geom_point(aes(x=reorder(Gene_Symbol, meanControl), y=-log10(meanControl)), size=1, colour=ggcol[3], position=position_jitter(width=.2, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, meanControl), y=-log10(meanAvrBs3)), size=1, colour=ggcol[1], position=position_jitter(width=.2, height=0)) +
            geom_point(aes(x=reorder(Gene_Symbol, meanControl), y=-log10(meanAvrBs4)), size=1, colour=ggcol[2], position=position_jitter(width=.2, height=0)) +
            ylab("Mean TPM") +
            xlab("Gene symbols") +
            plot_theme +
            theme(axis.text.x=element_text(angle=90, size=8, vjust=0.4)) +
            scale_y_continuous(breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
        grid.newpage()
        grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
        dev.off()

        table_EBEs_both = both_merged[, c(1, ncol(both_merged), ncol(both_merged)-2, ncol(both_merged)-1)]
        setDT(table_EBEs_both)[setDT(test_res_both), meanControl := i.meanControl, on=c("Gene_Symbol")]
        setDT(table_EBEs_both)[setDT(test_res_both), meanAvrBs3 := i.meanAvrBs3, on=c("Gene_Symbol")]
        setDT(table_EBEs_both)[setDT(test_res_both), meanAvrBs4 := i.meanAvrBs4, on=c("Gene_Symbol")]
        table_EBEs_both = cbind(table_EBEs_both, both_merged[, c(2, 3, 4, 5, 6, 7)])
        table_EBEs_both = table_EBEs_both[with(table_EBEs_both, order(-AvrBs3_logFC)), ]
        rownames(table_EBEs_both) = NULL
        write.table(table_EBEs_both , file=paste(maindir, "/", reg, "/DEG_illustrations_tables/DEG_both_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

        table_EBEs_both_pred = table_EBEs_both[grep(paste(overlap_both, collapse="|"), table_EBEs_both$Gene_Symbol), ]
        write.table(table_EBEs_both_pred , file=paste(maindir, "/", reg, "/DEG_illustrations_tables/DEG_both_EBE_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

        table_EBEs_both_nopred = table_EBEs_both[!grep(paste(overlap_both, collapse="|"), table_EBEs_both$Gene_Symbol), ]
        write.table(table_EBEs_both_nopred , file=paste(maindir, "/", reg, "/DEG_illustrations_tables/DEG_both_noEBE_table.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)


        ################################################################


        data_edge_TPM_wide = data_edge_TPM
        data_edge_TPM_wide$Length = NULL
        data_edge_TPM_long = reshape(data_edge_TPM_wide, varying = 2:ncol(data_edge_TPM_wide), sep = "_", direction = 'long')
        data_edge_TPM_long = melt(data_edge_TPM_long)

        options(ucscChromosomeNames=FALSE)
        grtrack <- GeneRegionTrack(genome_annot, name="Annotation")
        gtrack <- GenomeAxisTrack()

        main_bams = paste(maindir, "/", "mappings", sep="")
        #bam_sets = c("AvrBs3_merged.bam", "AvrBs4_merged.bam", "empty_merged.bam")
        treatments = c("AvrBs3", "AvrBs4", "empty")

        inputs = c("AvrBs3", "AvrBs4")
        for(Avr in inputs){

            if(Avr == "AvrBs3"){
                #myset = names(AvrBs3_EBEs_orig)[2:length(names(AvrBs3_EBEs_orig))]
                bam_sets = c("AvrBs3_merged.bam", "empty_merged.bam")
                treatments = c("AvrBs3", "empty")
                hitsEBE = table_EBEs_AvrBs3_pred
                hitsEBE_symbols = as.vector(hitsEBE[, hitsEBE$Gene_Symbol])

                # keep = c()
                # for(count in seq(1, length(AvrBs3_EBEs_orig), by=1)){
                #     EBE_DNAStringSet = AvrBs3_EBEs_orig[count]
                #     name.dist = names(EBE_DNAStringSet)
                #     name.dist_split = strsplit(name.dist, "\\.")
                #     if(name.dist_split[[1]][1] %in% hitsEBE_symbols){
                #         keep = c(keep, count)
                #     }
                # }
                hitsEBE_EBEs = AvrBs3_EBEs_orig
                hitsEBE_ids = names(AvrBs3_EBEs_orig)
                # hitsEBE_EBEs = append(TALE_AvrBs3, AvrBs3_EBEs_orig[keep])

                hitsEBE_table = table_EBEs_AvrBs3_pred[grepl(paste(hitsEBE_symbols, collapse="|"), table_EBEs_AvrBs3_pred$Gene_Symbol),]
                pdf(paste(maindir, "/", reg, "/TSS/images/AvrBs3_", reg, "reg_predEBEs.pdf", sep=""), 15, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
                p1 = ggplot(hitsEBE_table) +
                    geom_bar(aes(x=reorder(Gene_Symbol, -meanAvrBs3), y=(AvrBs3_logFC)), stat="identity", colour="black", fill=ggcol[1]) +
                    ylab("Fold change") +
                    plot_theme +
                    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
                    scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
                p2 = ggplot(hitsEBE_table) +
                    geom_point(aes(x=reorder(Gene_Symbol, -meanAvrBs3), y=-log10(meanControl)), size=3, colour=ggcol[3]) +
                    geom_point(aes(x=reorder(Gene_Symbol, -meanAvrBs3), y=-log10(meanAvrBs3)), size=3, colour=ggcol[1]) +
                    ylab("Mean TPM") +
                    xlab("Gene symbols") +
                    plot_theme +
                    theme(axis.text.x=element_text(angle=90, size=9, vjust=0.4)) +
                    scale_y_continuous(limits=c(-3, 2), breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
                grid.newpage()
                grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
                dev.off()


            } else{
                #myset = names(AvrBs4_EBEs_orig)[2:length(names(AvrBs4_EBEs_orig))]
                bam_sets = c("AvrBs4_merged.bam", "empty_merged.bam")
                treatments = c("AvrBs4", "empty")
                hitsEBE = table_EBEs_AvrBs4_pred
                hitsEBE_symbols = as.vector(hitsEBE[, hitsEBE$Gene_Symbol])

                # keep = c()
                # for(count in seq(1, length(AvrBs4_EBEs_orig), by=1)){
                #     EBE_DNAStringSet = AvrBs4_EBEs_orig[count]
                #     name.dist = names(EBE_DNAStringSet)
                #     name.dist_split = strsplit(name.dist, "\\.")
                #     if(name.dist_split[[1]][1] %in% hitsEBE_symbols){
                #         keep = c(keep, count)
                #     }
                # }
                hitsEBE_EBEs = AvrBs4_EBEs_orig
                hitsEBE_ids = names(AvrBs4_EBEs_orig)

                hitsEBE_table = table_EBEs_AvrBs4_pred[grepl(paste(hitsEBE_symbols, collapse="|"), table_EBEs_AvrBs4_pred$Gene_Symbol),]
                pdf(paste(maindir, "/", reg, "/TSS/images/AvrBs4_", reg, "reg_predEBEs.pdf", sep=""), 7, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
                p1 = ggplot(hitsEBE_table) +
                    geom_bar(aes(x=reorder(Gene_Symbol, -meanAvrBs4), y=(AvrBs4_logFC)), stat="identity", colour="black", fill=ggcol[2]) +
                    ylab("Fold change") +
                    plot_theme +
                    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
                    scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
                p2 = ggplot(hitsEBE_table) +
                    geom_point(aes(x=reorder(Gene_Symbol, -meanAvrBs4), y=-log10(meanControl)), size=3, colour=ggcol[3]) +
                    geom_point(aes(x=reorder(Gene_Symbol, -meanAvrBs4), y=-log10(meanAvrBs4)), size=3, colour=ggcol[2]) +
                    ylab("Mean TPM") +
                    xlab("Gene symbols") +
                    plot_theme +
                    theme(axis.text.x=element_text(angle=90, size=9, vjust=0.4)) +
                    scale_y_continuous(limits=c(-3, 2), breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
                grid.newpage()
                grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
                dev.off()
            }


            print("")
            print(Avr)

            # TSS statistics-related
            stats_dir = paste(maindir, "/", reg, "/TSS/statistics/", Avr, sep="")
            dir.create(stats_dir, showWarnings=FALSE)
            dir.create(paste(maindir, "/", reg, "/TSS/statistics/", Avr, "/tmp", sep=""), showWarnings=FALSE)
            dir.create(paste(maindir, "/", reg, "/TSS/images/", Avr, sep=""), showWarnings=FALSE)

            treatment_TPMs = data_edge_TPM[, c(1, as.vector(grep(Avr, names(data_edge_TPM))))]
            treatment_TPMs[] = lapply(treatment_TPMs, function(x) type.convert(as.character(x)))

            treatment_TPMs$mean = rowMeans(treatment_TPMs[,2:ncol(treatment_TPMs)])
            treatment_TPMs = treatment_TPMs[order(treatment_TPMs$mean), ]
            no_groups = 20
            group_size = round(nrow(treatment_TPMs)/no_groups, 0)

            diff = nrow(treatment_TPMs) - (group_size * no_groups)
            dist_diff = ceil(diff/no_groups)

            groups_col = list()
            for(n in 1:no_groups){
                group = rep(n, group_size + dist_diff)
                groups_col = append(groups_col, group)

                diff = diff - dist_diff
                if(diff == 0){
                    dist_diff = 0
                }
            }

            treatment_TPMs$group_TPM = as.vector(unlist(groups_col))

            treatment_TPMs_file = paste(maindir, "/", reg, "/TSS/", Avr, "_treatment_TPMs.tsv", sep="")
            write.table(treatment_TPMs, file=treatment_TPMs_file, row.names=FALSE, quote=FALSE, sep="\t")

            ggplot(data=treatment_TPMs, aes(x=factor(treatment_TPMs$group_TPM), y=treatment_TPMs$mean)) +
                geom_boxplot() +
                ylab("TPMs") +
                xlab("Group number") +
                scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels)) +
                plot_theme
            setwd(stats_dir)
            ggsave(paste(Avr, "_groups.pdf", sep=""), width=20, units="cm")
            ggsave(paste(Avr, "_groups.png", sep=""), width=20, units="cm", dpi=450)

            if(Avr == "AvrBs3"){
                library_size_group = "2_AvrBs3"
            } else{
                library_size_group = "3_AvrBs4"
            }

            library_sizes = as.data.frame(edgeR_dge$samples)
            library_size_control = sum(library_sizes[library_sizes["group"] == "1_control", ]["lib.size"] * library_sizes[library_sizes["group"] == "1_control", ]["norm.factors"])
            library_size_treatment = sum(library_sizes[library_sizes["group"] == library_size_group, ]["lib.size"] * library_sizes[library_sizes["group"] == library_size_group, ]["norm.factors"])
            ######################


            error_file = paste(maindir, "/", reg, "/TSS/images/", Avr, "_errors.txt", sep="")
            if(file.exists(error_file)){
                file.remove(error_file)
            }

            count_EBEs = 1
            selected_genes = list()
            selected_genes_EBEs = list()
            indgene_list = list()
            n_list = list()
            z_list = list()
            pnorm_list = list()
            pwilcox_list = list()

            for(indgene in hitsEBE_ids){
                print(Sys.time())
                print(indgene)
                percent = round((count_EBEs/length(hitsEBE_ids)) * 100, 1)
                count_EBEs = count_EBEs + 1
                print(paste("Overall progress: ", percent, " %", sep=""))

                indgene_split = strsplit(indgene, "\\.")
                indgene_name = indgene_split[[1]][1]

                indgene_EBEpos = as.numeric(indgene_split[[1]][2])

                bams_dists = c(0, 0)
                ii = 1

                pdf(paste(maindir, "/", reg, "/TSS/images/", Avr, "/", indgene, ".pdf", sep=""), 9, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")

                for(bam_name in bam_sets){
                    distance = "?"

                    treatment = treatments[ii]

                    bamfile_path = paste(main_bams, "/", bam_name, sep="")

                    altrack <- AlignmentsTrack(bamfile_path, name="Mapped reads")
                    displayPars(altrack) = list(min.height=0)

                    indgene_gff = gff_gene[gff_gene$gene == indgene_name, ]
                    indgene_chrom = as.character(seqnames(indgene_gff))
                    indgene_gff_start = as.numeric(start(indgene_gff))
                    indgene_gff_end = as.numeric(end(indgene_gff))
                    indgene_gff_strand = as.character(strand(indgene_gff))

                    # indgene_gff = gff_cds[gff_cds$gene == indgene_name, ]
                    # indgene_gff_cds = unique(indgene_gff$ID)
                    #
                    # indgene_gff_useCDS = indgene_gff_cds[1]
                    #
                    # if(is.na(indgene_gff_useCDS)){
                    #     indgene_gff = gff_gene[gff_gene$gene == indgene_name, ]
                    #     indgene_chrom = as.character(seqnames(indgene_gff))
                    #     indgene_gff_start = as.numeric(start(indgene_gff))
                    #     indgene_gff_end = as.numeric(end(indgene_gff))
                    #     indgene_gff_strand = as.character(strand(indgene_gff))
                    # } else{
                    #     indgene_gff = gff_cds[gff_cds$ID == indgene_gff_useCDS, ]
                    #     indgene_gff_useCDS_starts = start(indgene_gff)
                    #     indgene_gff_useCDS_ends = end(indgene_gff)
                    #     indgene_chrom = as.character(seqnames(indgene_gff)[1])
                    #     indgene_gff_start = as.numeric(min(indgene_gff_useCDS_starts))
                    #     indgene_gff_end = as.numeric(max(indgene_gff_useCDS_ends))
                    #     indgene_gff_strand = as.character(strand(indgene_gff)[1])
                    # }

                    if(indgene_gff_strand == "+"){
                        indgene_start = indgene_gff_start + indgene_EBEpos
                        indgene_end = indgene_start + 300
                        highlight = HighlightTrack(trackList=list(grtrack), start=indgene_start, width=18, chromosome=indgene_chrom)
                    } else if(indgene_gff_strand == "-"){
                        indgene_end = indgene_gff_end - indgene_EBEpos
                        indgene_start = indgene_end - 300
                        highlight = HighlightTrack(trackList=list(grtrack), start=indgene_end-18, width=18, chromosome=indgene_chrom)
                    }

                    gene_reads = GRanges(seqnames=indgene_chrom, IRanges(indgene_start, indgene_end))
                    what_info = c("rname", "strand", "pos", "qwidth")
                    param = ScanBamParam(which=gene_reads, what=what_info, simpleCigar=TRUE)
                    bamFile = BamFile(bamfile_path, index=paste(bamfile_path, ".", "bai", sep=""))

                    bam = scanBam(bamFile, param=param)
                    reads_pos = as.character(unlist(bam[[1]][3]))
                    reads_length = as.character(unlist(bam[[1]][4]))

                    size_reads = as.numeric(length(reads_pos))

                    if(size_reads != 0){

                        if(indgene_gff_strand == "+"){
                            min_pos = which.min(reads_pos)
                            distance = -1 * (as.numeric(reads_pos[min_pos]) - indgene_start - 18)

                        } else if(indgene_gff_strand == "-"){
                            max_pos = which.max(reads_pos)
                            distance = -1 * ((indgene_end - 18) - (as.numeric(reads_pos[max_pos]) + as.numeric(reads_length[max_pos])))
                        }
                    }

                    if(distance > 0){
                        distance = "x"
                    }


                    if(distance > -150 | distance == "x" | ii == 2 | reg == "down"){
                        bams_dists[ii] = 1

                        selected_genes = append(selected_genes, indgene_name)
                        selected_genes_EBEs = append(selected_genes_EBEs, indgene)

                    }

                    if(bams_dists[1] == 1){

                        if(!(grepl("empty", bam_name)) & reg != "down"){
                            testGene = NA
                            z_score = "x"
                            norm_p = "x"
                            wilcox_p = "x"

                            TSS_file = paste(stats_dir, "/", indgene, ".tsv", sep="")

                            TSS_cmd = paste("/home/sven/Documents/Python/Projects/enrichTSS/enrichTSS_multi.py -@ 6 -n 18 -w 55 -t ", indgene_name, " -e ", indgene_EBEpos,
                                            " -x ", treatment_TPMs_file, " -g ", mod_annot, " -m /home/sven/Documents/Uni/Masterthesis/Analyses/2017_ECW-123_Xcv_AvrBs3_AvrBs4/Bioconductor/20170507093925/mappings -s ", paste(bam_sets, collapse=" "), " -o ", TSS_file, " -c ", Avr, " --libsize ", library_size_treatment, " ", library_size_control, " --is_pe", sep="")
                            system(TSS_cmd)

                            comp_group = read.table(TSS_file, header=TRUE)
                            comp_group$cpm_enrichment_scaled = scale(log10(comp_group$cpm_enrichment))
                            testGene = comp_group[comp_group$Gene_IDs == indgene_name, ]

                            if(length(testGene$Gene_IDs) > 0 & nrow(comp_group) >= 30){
                                comp_mean = mean(comp_group$cpm_enrichment_scaled)
                                comp_sd = sd(comp_group$cpm_enrichment_scaled)

                                z_score = round((testGene$cpm_enrichment_scaled - comp_mean) / comp_sd, 3)

                                norm_p = round(pnorm(-abs(z_score)), 6)

                                wilcox_p = wilcox.test(as.vector(comp_group$cpm_enrichment_scaled), mu=as.numeric(testGene$cpm_enrichment_scaled), alternative="less")
                                wilcox_p = round(wilcox_p$p.value, 6)
                            }

                            ggplot(data=comp_group, aes(x=comp_group$cpm_enrichment_scaled)) +
                                geom_histogram(bins=50) +
                                plot_theme +
                                xlab("log10(enrichment factor)")
                            setwd(stats_dir)
                            ggsave(paste(indgene, ".pdf", sep=""), width=15, units="cm")

                            write.table(comp_group, file=paste(stats_dir, "/", indgene, "_n-", nrow(comp_group), "_z-", z_score, "_normp-", norm_p, "_wilcoxp-", wilcox_p, ".tsv", sep=""), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

                            file.remove(TSS_file)

                            indgene_list = append(indgene_list, indgene)
                            n_list = append(n_list, nrow(comp_group))
                            z_list = append(z_list, z_score)
                            pnorm_list = append(pnorm_list, norm_p)
                            pwilcox_list = append(pwilcox_list, wilcox_p)
                        }

                        ploterror = tryCatch({
                            plotTracks(list(gtrack, highlight, altrack), from=indgene_start, to=indgene_end, chromosome=indgene_chrom, littleTicks=TRUE, exponent=0, main=paste(treatment, "; ", indgene_name, "; dist. EBE to 5' end of gene: ", indgene_EBEpos, ", dist. EBE to 1st transcript: ", distance, sep=""))
                        }, error = function(err){
                            print("-- Plot error")
                            write.table(indgene, file=error_file, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
                        })

                    }

                    ii = ii + 1
                }
                dev.off()
                gc()
            }

            pdfs = list.files(path=paste(maindir, "/", reg, "/TSS/images/", Avr, sep=""), full.names=TRUE)
            for(file in pdfs){
                info = file.info(file)
                if(info$size <= 812){
                    file.remove(file)
                }
            }

            if(reg != "down"){
                indgene_df = data.frame(unlist(indgene_list), unlist(n_list), unlist(z_list), unlist(pnorm_list), unlist(pwilcox_list))
                names(indgene_df) = c("gene_EBE", "n", "z-score", "p_norm", "p_wilcox")
                write.table(indgene_df, file=paste(stats_dir, "/", Avr, "_statistics_results.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
            }

            subset_EBEs = hitsEBE_EBEs[unique(grep(paste(selected_genes_EBEs, collapse="|"), names(hitsEBE_EBEs)))]
            writeXStringSet(subset_EBEs, filepath=paste(maindir, "/", reg, "/TSS/images/", Avr, "/", Avr, "_EBEs.fa", sep=""), format="fasta")


            if(Avr == "AvrBs3"){
                hitsEBE_table = table_EBEs_AvrBs3_pred[grepl(paste(selected_genes, collapse="|"), table_EBEs_AvrBs3_pred$Gene_Symbol),]
                pdf(paste(maindir, "/", reg, "/TSS/images/AvrBs3_", reg, "reg_predEBEs_selection.pdf", sep=""), 11, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
                p1 = ggplot(hitsEBE_table) +
                    geom_bar(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=(AvrBs3_logFC)), stat="identity", colour="black", fill=ggcol[1]) +
                    ylab("Fold change") +
                    plot_theme +
                    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
                    scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
                p2 = ggplot(hitsEBE_table) +
                    geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanControl)), size=2, colour=ggcol[3]) +
                    geom_point(aes(x=reorder(Gene_Symbol, -AvrBs3_logFC), y=-log10(meanAvrBs3)), size=2, colour=ggcol[1]) +
                    ylab("Mean TPM") +
                    xlab("Gene symbols") +
                    plot_theme +
                    theme(axis.text.x=element_text(angle=90, size=9, vjust=0.4)) +
                    scale_y_continuous(limits=c(-3, 2), breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
                grid.newpage()
                grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
                dev.off()

            } else if(Avr == "AvrBs4"){
                hitsEBE_table = table_EBEs_AvrBs4_pred[grepl(paste(selected_genes, collapse="|"), table_EBEs_AvrBs4_pred$Gene_Symbol),]
                pdf(paste(maindir, "/", reg, "/TSS/images/AvrBs4_", reg, "reg_predEBEs_selection.pdf", sep=""), 7, 5, pointsize=0.5, pagecentre=TRUE, colormodel="rgb")
                p1 = ggplot(hitsEBE_table) +
                    geom_bar(aes(x=reorder(Gene_Symbol, -AvrBs4_logFC), y=(AvrBs4_logFC)), stat="identity", colour="black", fill=ggcol[2]) +
                    ylab("Fold change") +
                    plot_theme +
                    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
                    scale_y_continuous(trans="log10", breaks=breaks, labels=parse(text=labels))
                p2 = ggplot(hitsEBE_table) +
                    geom_point(aes(x=reorder(Gene_Symbol, -AvrBs4_logFC), y=-log10(meanControl)), size=2, colour=ggcol[3]) +
                    geom_point(aes(x=reorder(Gene_Symbol, -AvrBs4_logFC), y=-log10(meanAvrBs4)), size=2, colour=ggcol[2]) +
                    ylab("Mean TPM") +
                    xlab("Gene symbols") +
                    plot_theme +
                    theme(axis.text.x=element_text(angle=90, size=9, vjust=0.4)) +
                    scale_y_continuous(limits=c(-3, 2), breaks=c(3, 2, 1, 0, -1, -2), labels=c(0.001, 0.01, 0.1, 1, 10, 100))
                grid.newpage()
                grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
                dev.off()
            }
        }
    }
}



stats = read.table(file="/home/sven/Documents/Uni/Masterthesis/Analyses/2017_ECW-123_Xcv_AvrBs3_AvrBs4/Bioconductor/20170720015916/up/TSS/statistics/AvrBs4/AvrBs4_statistics_results.tsv", header=TRUE)

test = test_res_AvrBs4
test = test[test$AvrBs4 == 1, ]
test$Gene_Symbol = as.character(test$Gene_Symbol)
test = test[order(-test$AvrBs4_logFC), ]
stats = stats[order(stats$p_norm), ]
stats_genes = as.character(stats$gene_EBE)

for(ii in 1:length(stats_genes)){
    x = strsplit(stats_genes[ii], "\\.")
    stats_genes[ii] = x[[1]][1]
}

stats_genes = unique(stats_genes)
stats_genes_df = as.data.frame(stats_genes)
stats_genes_df$num = as.numeric(seq(1:nrow(stats_genes_df)))

test_genes_red = c()
test_genes_num = c()
for(gene in test$Gene_Symbol){
    if(gene %in% stats_genes){
        test_genes_red = c(test_genes_red, gene)

        num = stats_genes_df[stats_genes_df$stats_genes == gene, ]
        test_genes_num = c(test_genes_num, as.numeric(num$num))
    }
}

test_genes_df = as.data.frame(test_genes_red)
test_genes_df$num = test_genes_num

print(cor(stats_genes_df$num, test_genes_df$num, method="pearson"))


tpms = read.table(file="/home/sven/Documents/Uni/Masterthesis/Analyses/2017_ECW-123_Xcv_AvrBs3_AvrBs4/Bioconductor/20170720015916/up/TSS/AvrBs4_treatment_TPMs.tsv", header=TRUE)
tpms = tpms[order(-tpms$mean), ]
tpms_genes_red = c()
tpms_genes_num = c()
for(gene in tpms$Gene_IDs){
    if(gene %in% stats_genes){
        tpms_genes_red = c(tpms_genes_red, gene)

        num = stats_genes_df[stats_genes_df$stats_genes == gene, ]
        tpms_genes_num = c(tpms_genes_num, as.numeric(num$num))
    }
}
tpms_genes_df = as.data.frame(tpms_genes_red)
tpms_genes_df$num = tpms_genes_num

print(cor(stats_genes_df$num, tpms_genes_df$num, method="pearson"))



# AvrBs3_all_WTallDeriv = AvrBs3_merged[(rowSums(is.na(AvrBs3_merged[2:ncol(AvrBs3_merged)-1])) == 0),]
# AvrBs3_TALVEZ_WTallDeriv = AvrBs3_merged[(rowSums(is.na(AvrBs3_merged[2:8])) == 0),]
# AvrBs3_TALgetter_WTallDeriv = AvrBs3_merged[(rowSums(is.na(AvrBs3_merged[9:15])) == 0),]
# AvrBs3_TALENT_WTallDeriv = AvrBs3_merged[(rowSums(is.na(AvrBs3_merged[16:22])) == 0),]
#
# AvrBs3_merged = AvrBs3_merged_full
# pred_pos1 = c(3, 10, 17)
# pred_pos2 = c(8, 15, 22)
# tale = c("TALVEZ_AvrBs3", "TALgetter_AvrBs3", "TALENT_AvrBs3")
# target_num = c(0, 1, 2, 3, 4, 5)
# for(ii in 1:3){
#     for(num in target_num){
#     Deriv = AvrBs3_merged[(rowSums(is.na(AvrBs3_merged[pred_pos1[ii]:pred_pos2[ii]])) == num) & (is.na(AvrBs3_merged[tale[ii]]) == 1),]
#     Deriv = Deriv[, c(1, grep("AvrBs3$", names(Deriv), perl=T), grep("2xrep17$", names(Deriv), perl=T), grep("5\\.NS$", names(Deriv), perl=T), grep("drep9$", names(Deriv), perl=T), grep("drep15\\.16$", names(Deriv), perl=T), grep("7NS\\.drep11\\.13$", names(Deriv), perl=T), grep("drep16$", names(Deriv), perl=T), 23, 24) ]
#
#     if(nrow(Deriv) != 0){
#         write.table(Deriv , file=paste(maindir, "/TALE_EBEs/results/DerivWT/2017_EBEs_", tale[ii], "_", num, "-Deriv.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
#     }
#
#     WTDeriv = AvrBs3_merged[(rowSums(is.na(AvrBs3_merged[pred_pos1[ii]:pred_pos2[ii]])) == num) & (is.na(AvrBs3_merged[tale[ii]]) == 0),]
#     WTDeriv = Deriv[, c(1, grep("AvrBs3$", names(WTDeriv), perl=T), grep("2xrep17$", names(WTDeriv), perl=T), grep("5\\.NS$", names(WTDeriv), perl=T), grep("drep9$", names(WTDeriv), perl=T), grep("drep15\\.16$", names(WTDeriv), perl=T), grep("7NS\\.drep11\\.13$", names(WTDeriv), perl=T), grep("drep16$", names(WTDeriv), perl=T), 23, 24) ]
#     if(nrow(WTDeriv) != 0){
#         write.table(WTDeriv , file=paste(maindir, "/TALE_EBEs/results/DerivWT/2017_EBEs_", tale[ii], "_", num, "-Deriv_WT.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
#     }
#     }
# }
#
#
#
# AvrBs4_all_WTallDeriv = AvrBs4_merged[(rowSums(is.na(AvrBs4_merged[2:ncol(AvrBs4_merged)-1])) == 0),]
# AvrBs4_TALVEZ_WTallDeriv = AvrBs4_merged[(rowSums(is.na(AvrBs4_merged[grep("TALVEZ", names(AvrBs4_merged))])) == 0),]
# AvrBs4_TALgetter_WTallDeriv = AvrBs4_merged[(rowSums(is.na(AvrBs4_merged[grep("TALgetter", names(AvrBs4_merged))])) == 0),]
# AvrBs4_TALENT_WTallDeriv = AvrBs4_merged[(rowSums(is.na(AvrBs4_merged[grep("TALENT", names(AvrBs4_merged))])) == 0),]
#
# AvrBs4_merged = AvrBs4_merged_full
# pred_pos1 = c(2, 6)
# pred_pos2 = c(5, 9)
# tale = c("TALVEZ_AvrBs4", "TALgetter_AvrBs4")
# target_num = c(0, 1, 2, 3, 4)
# for(ii in 1:2){
#     for(num in target_num){
#         Deriv = AvrBs4_merged[(rowSums(is.na(AvrBs4_merged[pred_pos1[ii]:pred_pos2[ii]])) == num) & (is.na(AvrBs4_merged[tale[ii]]) == 1),]
#         Deriv = Deriv[, c(1, grep("AvrBs4$", names(Deriv), perl=T), grep("5\\.NI$", names(Deriv), perl=T), grep("drep15$", names(Deriv), perl=T), grep("drep13$", names(Deriv), perl=T),  13, 14) ]
#         if(nrow(Deriv) != 0){
#             write.table(Deriv , file=paste(maindir, "/TALE_EBEs/results/DerivWT/2017_EBEs_", tale[ii], "_", num, "-Deriv.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
#         }
#
#         WTDeriv = AvrBs4_merged[(rowSums(is.na(AvrBs4_merged[pred_pos1[ii]:pred_pos2[ii]])) == num) & (is.na(AvrBs4_merged[tale[ii]]) == 0),]
#         WTDeriv = WTDeriv[, c(1, grep("AvrBs4$", names(WTDeriv), perl=T), grep("5\\.NI$", names(WTDeriv), perl=T), grep("drep15$", names(WTDeriv), perl=T), grep("drep13$", names(WTDeriv), perl=T), 13, 14) ]
#         if(nrow(WTDeriv) != 0){
#             write.table(WTDeriv , file=paste(maindir, "/TALE_EBEs/results/DerivWT/2017_EBEs_", tale[ii], "_", num, "-Deriv_WT.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
#         }
#     }
# }
#
# # Other Analyses
# reg_AvrBs3_up_AvrBs4_down = test_res[test_res$AvrBs3 == 1 & test_res$AvrBs4 == -1, ]
# reg_AvrBs3_down_AvrBs4_up =  test_res[test_res$AvrBs3 == -1 & test_res$AvrBs4 == 1, ]
# reg_AvrBs3_some_AvrBs4_none = test_res[(test_res$AvrBs3 == -1 | test_res$AvrBs3 == 1) & test_res$AvrBs4 == 0, ]
# reg_both_down =  test_res[test_res$AvrBs3 == -1 & test_res$AvrBs4 == -1, ]

