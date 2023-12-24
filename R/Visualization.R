####Visualization########

#' PatternHeatmap: Heatmap to show the pattern specific temporal genes. Please first install Seurat, ggplot2, ComplexHeatmap and circlize packages.
#' 
#' @param obj The results of TDEseq analysis
#' @param features Genes to be shown in heatmap
#' @param features.show Genes to be annotated in heatmap
#' @param features.num Number of genes to be shown for each patterns 
#' @param cols Color of the heatmap.
#' @author Yue Fan, Shiquan Sun
#' @export
PatternHeatmap<-function(obj,features=NULL,features.show=NULL,features.num=50,cols=c("navy", "white", "firebrick3"))
{
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))
seu<-CreateSeuratObject(obj@NormalizeData)
seu@assays$RNA@data<-obj@NormalizeData
seu<-ScaleData(seu)
mat <- GetAssayData(seu,slot = 'scale.data')
metadata <- obj@Metadata
res_dat=obj@dfTDEseqResults
res_dat=res_dat[order(res_dat$pval),]
feature_annot=c()
if(is.null(features))
{
features_plot=c()
for(pattern in c("Growth","Recession","Peak","Trough"))
{
idx=which(res_dat$pattern==pattern)
if(length(idx)>=features.num)
{
features_plot=c(features_plot,res_dat$gene[idx[1:features.num]])
feature_annot=c(feature_annot,rep(pattern,features.num))
}else{
features_plot=c(features_plot,res_dat$gene[idx])
feature_annot=c(feature_annot,rep(pattern,length(idx)))
}
}
}else{
features_plot=c()
idx=match(features,res_dat$gene)
if(any(is.na(idx)))
{
idx=idx[-which(is.na(idx))]
}
features=features[idx]
subres_dat=res_dat[idx,]
for(pattern in c("Growth","Recession","Peak","Trough"))
{
idx=which(subres_dat$pattern==pattern)
features_plot=c(features_plot,res_dat$gene[idx])
feature_annot=c(feature_annot,rep(pattern,length(idx)))
}

}

group_info <- metadata$Time_Origin
col_fun = colorRamp2(c(-2, 0, 2),cols)

mat=mat[features_plot,]

if(!is.null(features.show))
{
gene_pos <- match(features.show,rownames(mat))
row_anno <- rowAnnotation(gene=anno_mark(at=gene_pos,labels = features.show,
                                         labels_gp=gpar(fontface = 3)))
										 


f1=Heatmap(mat,
           name = 'Expression',
           col = col_fun,
           cluster_rows = FALSE,                 
           cluster_columns = F,                  
           show_column_dend=F,                    
           show_row_dend = F,                     
           show_row_names = FALSE,                   
           show_column_names = F,                 
           column_split = group_info,      
           row_split = feature_annot,		   
           right_annotation = row_anno,                  
           border_gp = gpar(col = "black", lty = 2) 
           )	
}else{

f1=Heatmap(mat,
           name = 'Expression',
           col = col_fun,
           cluster_rows = FALSE,                 
           cluster_columns = F,                  
           show_column_dend=F,                    
           show_row_dend = F,                     
           show_row_names = FALSE,                   
           show_column_names = F,                 
           column_split = group_info,      
           row_split = feature_annot,		                 
           border_gp = gpar(col = "black", lty = 2) 
           )	


}		   
return(f1)
}

#' PatternHeatmap: Feature plot to show the pattern specific temporal genes. Please first install Seurat, ggplot2 and tidydr packages.
#' 
#' @param obj The results of TDEseq analysis
#' @param features Genes to be shown in feature plot
#' @author Yue Fan, Shiquan Sun
#' @export
PatternFeature<-function(obj,feature)
{
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tidydr"))
#Set background and format
tsne_theme <- theme( 
  axis.line=element_blank(), 
  axis.text.x=element_blank(), 
  axis.text.y=element_blank(), 
  axis.ticks=element_blank(), 
  axis.title.x=element_blank(), 
  axis.title.y=element_blank(), 
  panel.background=element_blank(), 
  panel.border=element_blank(), 
  panel.grid.major=element_blank(), 
  panel.grid.minor=element_blank()) +
  theme_dr(xlength = 0.15, ylength = 0.15,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed"))+
  theme(panel.grid = element_blank())+theme(text = element_text(size=10))+
  theme(axis.line.x=element_line(size=1))+
  theme(axis.line.y=element_line(size=1))

p <- FeaturePlot(
  seurat,
  features=feature,
  cols = c('#690E61', '#28748D','#F6E529'),
  max.cutoff='q98',label.size = 20)
p=p+ tsne_theme
return(p)
}

#' PatternLine: Line plot to show the pattern specific temporal genes. Please first install ggplot2.
#' 
#' @param seuobj A seurat object with UMAP embeddings been calculated
#' @param features Genes to be shown in feature plot
#' @author Yue Fan, Shiquan Sun
#' @export
PatternLine<-function (obj, feature.show = NULL, cols = NULL) 
{
    if(is.null(cols))
	{
	cols=c('#E50C7D',"#E34627","#A22066","#A474A4","#2D8573","#E1DE15","#C16728","#2578BE","#738DC8","#C0C0C0", '#7d8d8e','#2a24d0','#a2292e','#274382','#838d36')
	}
    suppressPackageStartupMessages(library("ggplot2"))
	feature.show=feature.show[feature.show%in%rownames(obj@ModelFits)]
	if(length(feature.show)==0)
	{
	stop(paste0("Genes to be shown are not specified!!"))
	}
    dat = obj@ModelFits[feature.show, ]
    time = sort(unique(obj@Metadata$Time_Origin))
	N=length(feature.show)
	if(N!=1)
	{
	data=lapply(1:N,function(x){dat=data.frame(stage=time,expr=dat[x,],feature=feature.show[x])
	return(dat)
	})
	data=do.call(rbind,data)
	}else{
	data=data.frame(stage=time,expr=dat,feature=feature.show)
	}
    p <- ggplot(data = data, aes(x = stage, y = expr, color = feature)) + 
        stat_smooth(method = "loess", aes(col = feature), se = FALSE, size = 3) + 
		theme_classic() + 
		xlab("Stage") + 
		ylab("log(count+1)")+
		scale_colour_manual(values=cols)+
		theme(axis.title=element_text(size=rel(2),face="bold"),
              axis.text=element_text(size=rel(2),face="bold"),
              legend.title=element_blank(),
              legend.text=element_text(size=rel(2),face="bold"))
    return(p)
}