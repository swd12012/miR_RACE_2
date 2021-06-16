library(DESeq2)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
library(EnhancedVolcano)
library( "gplots" )
library( "RColorBrewer" )
library( "genefilter" )

#Read in sample coding table
sampleInfo <- read.table('design.txt', header=T)

#Read in count data
countdata <- read.table('data/out/counts.txt', header=T, row.names=1)

#Remove extraneous rows in countdata
countdata <- countdata[,6:10]

#rename and reorder column names
sampleNames <- c('P68', 'P69', 'P70', 'P71', 'P72')
names(countdata) <- sampleNames
countdata <- countdata[,c(4,5,1:3)]

#exclude P69
countdata2 <- countdata[,c(1:3,5)]

#Read in annotation file
#annotations <- read.table('mouseENSEMBLID_genename.txt',header=T,sep=',')

#annotated_countdata <- merge(countdata,annotations, by.x='row.names',by.y='Gene.stable.ID')

#Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=sampleInfo, design=~Condition)

dds <- DESeq(dds)
res <- results(dds)
res

write.csv(as.data.frame(res), file='data/out2/results.csv')

#Exploratory data analysis
#plotMA( res, ylim = c(-1, 1) )
#plotDispEsts( dds )
hist( res$pvalue, breaks=20, col="grey" )
###  throw out lowly expressed genes?? ... I leave this as an exercise
###  add external annotation to "gene ids"
# log transform
rld = rlog( dds )
## this is where you could just extract the log transformed normalized data for each sample ID, and then roll your own analysis
head( assay(rld) )
mydata = assay(rld)

#sampleDists = dist( t( assay(rld) ) )
# heat map
#sampleDistMatrix = as.matrix( sampleDists )
#rownames(sampleDistMatrix) = rld$TissueCode
#colnames(sampleDistMatrix) = NULL
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#heatmap.2( sampleDistMatrix, trace="none", col=colours)
# PCs
# wow you can sure tell tissue apart
#print( plotPCA( rld, intgroup = "TissueCode") )
# heat map with gene clustering
library( "genefilter" )
# these are the top genes (that tell tissue apart no doubt)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )

png('plots/heatmap.png', width=8, height=8, res=300, units='in')
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), cexRow=0.5, cexCol=0.5, margins=c(12,8),srtCol=45)
dev.off()
# volcano plot this is an exercise

#res <- res[order(res$padj),]

#deg=subset(res,padj<0.05)

#write.table(as.data.frame(deg), file='DEGs.txt')

#Lifted this function off the internet to turn ENSEMBL IDs to gene symbols
library(org.Mm.eg.db)
convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select( 
  db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

heatmap_matrix <- assay(rld)[topVarGenes,]
heatmap_matrix$gene_symbol <- convertIDs(row.names(heatmap_matrix), "ENSEMBL","SYMBOL", org.Mm.eg.db)

#Make volcano plot

png('plots/volcano.png', width=8, height=8, res=300, units='in')
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')
dev.off()
