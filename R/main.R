#' spatial_vector
#'
#' generating vectors for spatial relatioship speculator in one sample
#' 
#' @param object S4, a Seurat object
#' @param target character, a feature name (e.g., gene name in profile or a column name in meta.data). 
#' @param features vector, a collection of feature names, will be analyzed to explore their relationships with the target feature. 
#' @param objNM character, object name. Default, 'obj'. 
#' @param steps vector, a vector containing step sizes for simulating the displacement of spatial object. Default, seq(2,12,2). 
#' @param min.cutoffs numeric, a cutoff value on the feature to create a binary representation. If not specified (NULL), a mean value will be generated automatically.
#' @param threads numeric, number of threads for running, which is advised for large feature sets (≥50 features). Default, NULL. 
#' @param package.libpath character, the absolute file path to R library containing required packages (e.g., Seurat). This path is used when the threadsparameter is specified. Default, NULL. 
#' @param stat boolean, determines whether projected scores and magnitude values should be calculated for the features. Default, TRUE.  
#'
#' @return list, pool_df: a dataframe summarizing maximum and minimum delta Jaccard indexes for each step across all features. projected.score: a vector containing projected scores for the features. vector.len: a vector containing magnitude values for the features. target: target feature name. features: a collection of feature names.
#'
#' @export
spatial_vector = function(object,target,features,objNM='obj',steps=seq(2,12,2),min.cutoffs=NULL,threads=NULL,package.libpath=NULL,stat=TRUE){
  
  pool_list = lapply(1:length(steps),function(i){
    operator_steps = rep(steps[i],4)
    img_cord1 = suppressWarnings(spatial_adjust(object,slice_id='slice1',proj_name=names(object),rotate_mirror=c(0,'N'),feature=target,plot=F,verbose=F))
    map = c(rownames(as.matrix(Seurat::GetAssayData(object))),colnames(object@meta.data))
    
    pool = NULL
    if(is.null(threads)){
      pool = sapply(1:length(features),function(j){
        if(length(which(map %in% features[j]))>0){
          img_cord2 = suppressWarnings(spatial_adjust(object,slice_id='slice1',proj_name=names(object),rotate_mirror=c(0,'N'),feature=features[j],plot=F,verbose=F))
          tmp = suppressWarnings(spatial_cordstat(img_cord1,img_cord2,operator_steps=operator_steps,min.cutoffs=min.cutoffs,plot=F))
          range(tmp$mov_stat$diff_pct)
        }else{
          c(NA,NA)
        }
      })
    }else{
      library(parallel)
      library(Seurat)
      env = new.env()
      cl = makeCluster(threads)
      clusterExport(cl=cl,varlist=c("spatial_adjust","spatial_cordstat","package.libpath"),envir=env)
      if(!is.null(package.libpath)) clusterEvalQ(cl, .libPaths(package.libpath))
      # clusterEvalQ(cl,library(Seurat))
      pool = parSapply(cl,1:length(features),function(j){
        if(length(which(map %in% features[j]))>0){
          img_cord2 = suppressWarnings(spatial_adjust(object,slice_id='slice1',proj_name=names(object),rotate_mirror=c(0,'N'),feature=features[j],plot=F,verbose=F))
          tmp = suppressWarnings(spatial_cordstat(img_cord1,img_cord2,operator_steps=operator_steps,min.cutoffs=min.cutoffs,plot=F))
          range(tmp$mov_stat$diff_pct)
        }else{
          c(NA,NA)
        }
      })
      stopCluster(cl)
    }
    pool = as.data.frame(t(pool))
    
    colnames(pool) = c('min_djaccard','max_djaccard')
    rownames(pool) = paste0(objNM,'_',steps[i],'_',features)
    pool$feature = features
    pool$step = steps[i]
    pool$orig.ident = objNM
    pool
  })
  names(pool_list) = paste0('s',steps)
  pool_df = do.call(rbind,pool_list)
  
  projected.score = NULL
  vector.len = NULL
  if(stat){
    projected.score = spatial_vecProj(list(vectors=pool_df))
    vector.len = spatial_vecMagnitude(list(vectors=pool_df,features=features))
  }
  
  return(list(vectors=pool_df,projected.score=projected.score,vector.len=vector.len,target=target,features=features))
}

#' spatial_vectorX
#'
#' generating vectors for spatial relatioship speculator across various samples
#' 
#' @param objects list, a list containing a collection of Seurat objects. 
#' @param target character, a feature name (e.g., gene name in profile or a column name in meta.data). 
#' @param features vector, a collection of feature names, will be analyzed to explore their relationships with the target feature. 
#' @param steps vector, a vector containing step sizes for simulating the displacement of spatial object. Default, seq(2,12,2). 
#' @param min.cutoffs numeric, a cutoff value on the feature to create a binary representation. If not specified (NULL), a mean value will be generated automatically.
#' @param threads numeric, number of threads for running, which is advised for large feature sets (≥50 features). Default, NULL. 
#' @param package.libpath character, the absolute file path to R library containing required packages (e.g., Seurat). This path is used when the threadsparameter is specified. Default, NULL. 
#' @param output2RDS character, the absolute file path to outputting result. Default, NULL. 
#' 
#' @return list, pool_df: a dataframe summarizing maximum and minimum delta Jaccard indexes for each step across all features. pool_raw: a dataframe containing all delta Jaccard indexes for each step across all features. projected.score: a vector containing projected scores for the features. vector.len: a vector containing magnitude values for the features. target: target feature name. features: a collection of feature names. 
#'
#' @export
spatial_vectorX = function(objects,target,features,steps=seq(2,12,2),min.cutoffs=NULL,threads=NULL,package.libpath=NULL,output2RDS=NULL){
  
  pool_raw = c()
  for(i in 1:length(objects)){
    pool_raw = rbind(pool_raw,
                     spatial_vector(object=objects[[i]],
                                    target=target,
                                    features=features,
                                    objNM=names(objects)[i],
                                    steps=steps,
                                    min.cutoffs=min.cutoffs,
                                    threads=threads,
                                    stat=F,
                                    package.libpath=package.libpath)$vectors)
  }
  
  pool_df = as.data.frame(aggregate(.~feature+step,pool_raw[,1:4],mean))
  pool_df$orig.ident = 'obj'
  
  projected.score = spatial_vecProj(list(vectors=pool_df))
  vector.len = spatial_vecMagnitude(list(vectors=pool_df,features=features))
  
  final_obj = list(vectors=pool_df,pool_raw=pool_raw,projected.score=projected.score,vector.len=vector.len,target=target,features=features)
  if(!is.null(output2RDS)) saveRDS(final_obj,file=output2RDS)
  
  return(final_obj)
}

#' spatial_vecPlot
#' 
#' visualizing vectors for spatial relatioship speculator
#'
#' @param spaVecObj list, a spatial_vector or spatial_vectorX object. 
#' @param cex_point numeric, size of point. Default, 1. 
#' @param cex_text numeric, size of text. Default, 1. 
#' @param pct_ctf numeric, percentile cutoff for identifying features with extreme vector values, above which the feature name will be shown. Default, NULL. 
#' @param cols vector, a sequence of colors to use. If not specified (NULL), a default color gradient from 'steelblue' to 'red' will be generated automatically.
#' @param projected boolean, determines whether the projetced score should be visualized. Default, TRUE. 
#'
#' @return NULL
#' @export
spatial_vecPlot = function(spaVecObj,cex_point=1,cex_text=1,pct_ctf=NULL,cols=NULL,projected=T,...){
  
  cols_func = colorRampPalette(c('steelblue','orange','red'))
  
  spaVec.df = spaVecObj$vectors
  target = spaVecObj$target
  steps = sort(unique(spaVec.df$step))
  if(is.null(cols)) cols = cols_func(length(steps))
  rm(spaVecObj)
  
  init_step = as.data.frame(spaVec.df[which(spaVec.df$step == steps[1]),])
  final_step = as.data.frame(spaVec.df[which(spaVec.df$step == steps[length(steps)]),])
  
  xlim = range(spaVec.df$max_djaccard,na.rm=T)
  ylim = range(spaVec.df$min_djaccard,na.rm=T)
  position = apply(as.data.frame(final_step[,c('min_djaccard','max_djaccard')]),1,function(x){sum(x/2,na.rm=T)})
  xlim = c(min(c(xlim[1],position,0)),max(c(xlim[2],position,0)))
  ylim = c(min(c(ylim[1],position,0)),max(c(ylim[2],position,0)))
  
  op = par(no.readonly = TRUE)
  par(mar = c(5, 4, 4, 6) + 0.1)
  
  plot(init_step$max_djaccard,init_step$min_djaccard,xlim=xlim,ylim=ylim,cex=cex_point,col=cols[1],xlab='max djaccard',ylab='min djaccard',main=paste0('relation to ',target),...)
  if(projected) abline(a=0,b=1,lty=3)
  abline(v=0,lty=3); abline(h=0,lty=3)
  for(i in 2:length(steps)) points(spaVec.df[which(spaVec.df$step==steps[i]),'max_djaccard'],spaVec.df[which(spaVec.df$step==steps[i]),'min_djaccard'],cex=cex_point,col=cols[i],...)
  
  if(!is.null(pct_ctf)){
    x_max = quantile(final_step$max_djaccard[which(final_step$feature!=target)],c(pct_ctf,1-pct_ctf),na.rm = T)
    y_min = quantile(final_step$min_djaccard[which(final_step$feature!=target)],pct_ctf,na.rm = T)
    final_step$label = 0
    final_step$label[which(final_step$max_djaccard>x_max[2])] = 1
    final_step$label[which(final_step$max_djaccard<x_max[1])] = -1
    final_step$label[which(final_step$min_djaccard<y_min[1])] = -2
    feature_label= unique(final_step$feature[which(final_step$label!=0)])
    text(final_step$max_djaccard[which(final_step$feature %in% feature_label)],
         final_step$min_djaccard[which(final_step$feature %in% feature_label)],
         final_step$feature[which(final_step$feature %in% feature_label)],cex=cex_text)
  }else{
    text(final_step$max_djaccard,final_step$min_djaccard,final_step$feature,cex=cex_text)
  }
  
  for(i in 1:nrow(init_step)) lines(c(0,init_step$max_djaccard[i]),c(0,init_step$min_djaccard[i]),col=cols[1],...)
  
  for(i in 1:(length(steps)-1)){
    pre = as.data.frame(spaVec.df[which(spaVec.df$step == steps[i]),])
    pos = as.data.frame(spaVec.df[which(spaVec.df$step == steps[i+1]),])
    for(j in 1:nrow(pre)) lines(c(pre$max_djaccard[j],pos$max_djaccard[j]),c(pre$min_djaccard[j],pos$min_djaccard[j]),col=cols[i+1],...)
  } 
  
  if(projected) points(position,position,pch=3)
  
  par(xpd = TRUE)
  usr = par("usr")
  legend_width = 0.08 * (usr[2] - usr[1])  
  legend_x_left = usr[2] + 0.02 * (usr[2] - usr[1]) 
  legend_x_right = legend_x_left + legend_width
  legend_y_bottom = usr[3] 
  legend_y_top = usr[4] 
  
  for (i in 1:length(cols)) {
    rect(legend_x_left,
         legend_y_bottom + (i-1)/length(cols) * (legend_y_top - legend_y_bottom),
         legend_x_right,
         legend_y_bottom + i/length(cols) * (legend_y_top - legend_y_bottom),
         col = cols[i], border = NA)
  }
  
  z_range = range(steps)
  label_values = pretty(z_range) 
  label_values = label_values[which(label_values<=z_range[2])]
  label_positions = legend_y_bottom + (label_values - z_range[1]) / (z_range[2] - z_range[1]) * (legend_y_top - legend_y_bottom)
  axis(side = 4, at = label_positions, labels = label_values, las = 1, pos = legend_x_right, cex.axis = 0.8)
  
  par(mar = c(5, 4, 4, 2) + 0.1, xpd = FALSE)
  
}

#' spatial_magPlot
#' 
#' visualizing magnitude & projected score for spatial relatioship speculator
#'
#' @param spaVecObj list, a spatial_vector or spatial_vectorX object. 
#' @param cols vector, a sequence of colors to use. Default, c('green','orange'). 
#' @param cex_point numeric, size of point. Default, 1. 
#' @param cex_text numeric, size of text. Default, 1.  
#' @param cex_line numeric, size of line. Default, 1. 
#' @param lty numeric, line type. Default, 1. 
#'
#' @return NULL
#' @export
spatial_magPlot = function(spaVecObj,cols=c('green','orange'),cex_point=1,cex_text=1,cex_line=1,lty=1,...){
  
  spaVec.df = spaVecObj$vectors
  target = spaVecObj$target
  
  xvec = spaVecObj$vector.len
  yvec = spaVecObj$projected.score[names(xvec)]
  pos.shift = sd(xvec,na.rm=T)
  xlim = range(c(xvec,max(xvec,na.rm=T)+pos.shift),na.rm=T)
  
  plot(xvec,yvec,pch=3,xlim=xlim,
       xlab='magnitude',ylab='projected score',
       main=paste0('relation to ',target),
       cex=cex_point,...)
  abline(h=0,lty=3)
  text(xvec+pos.shift*0.5,yvec,label=names(xvec),cex=cex_text)
  for(i in 1:length(xvec)) lines(c(0,xvec[i]),rep(yvec[i],2),col=ifelse(yvec[i]>0,cols[1],cols[2]),cex=cex_line,lty=lty)
  
}

#' spatial_vecProj
#' 
#' calculating projected score
#'
#' @param spaVecObj list, a spatial_vector or spatial_vectorX object. 
#' @return vector, a vector containing projected scores for all specified features.
#' 
spatial_vecProj = function(spaVecObj){
  spaVec.df = spaVecObj$vectors
  steps = sort(unique(spaVec.df$step))
  projected.score = apply(as.data.frame(spaVec.df[which(spaVec.df$step == steps[length(steps)]),c('min_djaccard','max_djaccard')]),1,function(x){sum(x/2,na.rm=T)})
  names(projected.score) = spaVec.df[which(spaVec.df$step == steps[length(steps)]),'feature']
  return(projected.score)
}

#' spatial_vecMagnitude
#'
#' calculating magnitude
#'
#' @param spaVecObj list, a spatial_vector or spatial_vectorX object. 
#' @return vector, a vector containing magnitude values for specified features. 
#'
spatial_vecMagnitude = function(spaVecObj){
  spaVec.df = spaVecObj$vectors
  features = spaVecObj$features
  steps = sort(unique(spaVec.df$step))
  
  vector.len = NULL
  if(length(steps)>=2){
    vector.len = sapply(features,function(f){
      tmp = spaVec.df[which(spaVec.df$feature==f),c('min_djaccard','max_djaccard','step')]
      tmp1 = tmp[order(tmp$step),-3]
      tmp2 = rbind(c(0,0),tmp1[-nrow(tmp1),])
      sum(sqrt(rowSums((tmp1-tmp2)^2)),na.rm=T)
    })
  }else{
    vector.len = sapply(features,function(f){
      sqrt(rowSums(spaVec.df[which(spaVec.df$feature==f),c('min_djaccard','max_djaccard')]^2))
    }) 
  }
  
  return(vector.len)
}

#' spatial_adjust
#' 
#' creating a spatial object with digital matrix and showing spatial map with specified feature
#'
#' @param object S4, a Seurat object
#' @param slice_id character, slice id. Default, 'slice1'. 
#' @param proj_name character, a project name. Default, NULL. 
#' @param rotate_mirror vector, a vector of two elements controlling image transformation. The first element (rotate) is the rotation angle in degrees, valid values are 0, 90, 180, or 270. The second element (mirror) is a character specifying mirroring: 'Y' for yes, 'N' for no. Default is list(0, 'N'). 
#' @param feature character, a feature name (e.g., gene name in profile or a column name in meta.data). 
#' @param operator matrix, a user-defined 2x2 matrix for image transformation. This parameter is mutually exclusive with rotate_mirror. If both are provided, operator takes precedence. Default, NULL. 
#' @param cols vector, a sequence of colors to use. If not specified (NULL), a default color gradient from 'steelblue' to 'red' will be generated automatically. 
#' @param plot boolean, determines whether visulization should be generated. Default, TRUE. 
#' @param verbose boolean, determines whether running message should be outputting. Default, TRUE. 
#' @param pt_cex numeric, size of point. Default, 0.5. 
#'
#' @return list, img_cord: position values of the spots; vec_feature: feature values of the spots; vec_col: a sequence of colors; proj_name: project name; feature: feature name; operation: operator setting. 
#' @export
spatial_adjust = function(object,slice_id='slice1',proj_name=NULL,rotate_mirror=c(0,'N'),feature=NULL,operator=NULL,cols=NULL,plot=T, verbose=T, pt_cex=0.5){
  
  operator_cis0 = matrix(c(1,0,0,1),nrow=2)
  operator_cis90 = matrix(c(0,-1,1,0),nrow=2)
  operator_cis180 = matrix(c(-1,0,0,-1),nrow=2)
  operator_cis270 = matrix(c(0,1,-1,0),nrow=2)
  operator_mirror_N = matrix(c(1,0,0,1),nrow=2)
  operator_mirror_X = matrix(c(-1,0,0,1),nrow=2)
  operator_mirror_Y = matrix(c(1,0,0,-1),nrow=2)
  
  if(is.null(proj_name)) proj_name = slice_id
  if(length(which(rotate_mirror[1] %in% c(0,90,180,270)))==0) stop('[ERROR] rotate_mirror[1] should be 0, 90, 180, or 270. ')
  if(length(which(rotate_mirror[2] %in% c('N','X','Y')))==0) stop('[ERROR] rotate_mirror[2] should be "N", "X", or "Y". ')
  if(is.null(operator)) operator = eval(parse(text=paste0('operator_cis',rotate_mirror[1],' %*% operator_mirror_',rotate_mirror[2])))
  
  cols_func = colorRampPalette(c('steelblue','orange','red'))
  if(is.null(cols)) cols = cols_func(100)
  
  if(verbose) cat('[INFO] rotating',rotate_mirror[1],'degree, with',ifelse(rotate_mirror[2]=='N','no',ifelse(rotate_mirror[2]=='X','horizontal','vertical')),'mirror ..\n')
  raw_img = object@images[[slice_id]]@coordinates
  new_img = t(apply(raw_img[,c('row','col')],1,function(x,y){y %*% x},y=operator))
  rownames(new_img) = rownames(raw_img)
  colnames(new_img) = c('col_x','row_y')
  if(verbose) cat('[INFO] affecting',nrow(new_img),'spots ..\n')
  
  vec = NULL
  vec_col = NULL
  if(!is.null(feature)){
    if(verbose) cat('[INFO] mapping',feature,'into spots ..\n')
    mat = as.matrix(Seurat::GetAssayData(object))
    feature.idx = which(rownames(mat) == feature)
    vec = c()
    if(length(feature.idx)>0){
      vec = mat[feature.idx,]
      names(vec) = colnames(mat)
    } else if(any(colnames(object@meta.data)==feature)){
      vec = object@meta.data[[feature]]
      names(vec) = colnames(mat)
    } else{
      stop('[ERROR] invalid feature!')
    }
    rm(mat)
    
    if(length(cols)<2) stop('[ERROR] too few options for color mapping (length(cols) < 2)!')
    
    vec_factor = factor(round(vec,2),levels=sort(unique(round(vec,2))))
    vec_range = range(vec)
    
    if(verbose) cat('[INFO] the expression range of',feature,'is',vec_range,' ..\n')
    
    vec_col_idx = seq(vec_range[1],vec_range[2],length.out=length(cols)+1)
    vec_col = rep('#DDDDDD',length(vec))
    for(i in 1:length(vec_col_idx))  vec_col[which(vec>=vec_col_idx[i])] = cols[i]
  }
  
  if(plot){
    op = par(no.readonly = TRUE)
    par(mar = c(5, 4, 4, 6) + 0.1)
      
    plot(0,type='n',ylim=range(new_img[,2]),xlim=range(new_img[,1]),main=feature,xaxt='n',yaxt='n',xlab='',ylab='')
    for(i in 1:nrow(new_img)){
      if(!is.null(vec)) points(new_img[i,1],new_img[i,2],col=vec_col[i],pch=16,cex=pt_cex)
    }
      
    par(xpd = TRUE)
    usr = par("usr")
    legend_width = 0.08 * (usr[2] - usr[1])  
    legend_x_left = usr[2] + 0.02 * (usr[2] - usr[1]) 
    legend_x_right = legend_x_left + legend_width
    legend_y_bottom = usr[3] 
    legend_y_top = usr[4]
    
    for (i in 1:length(cols)) {
    rect(legend_x_left,
         legend_y_bottom + (i-1)/length(cols) * (legend_y_top - legend_y_bottom),
         legend_x_right,
         legend_y_bottom + i/length(cols) * (legend_y_top - legend_y_bottom),
         col = cols[i], border = NA)
    }
  
    z_range = range(vec)
    label_values = pretty(z_range) 
    label_values = label_values[which(label_values<=z_range[2])]
    label_positions = legend_y_bottom + (label_values - z_range[1]) / (z_range[2] - z_range[1]) * (legend_y_top - legend_y_bottom)
    axis(side = 4, at = label_positions, labels = label_values, las = 1, pos = legend_x_right, cex.axis = 0.8)
    par(mar = c(5, 4, 4, 2) + 0.1, xpd = FALSE)
    
  }
  
  return(list(
    img_cord = new_img,
    vec_feature = vec,
    vec_col = vec_col,
    proj_name = proj_name,
    feature = feature,
    operation = operator
  ))
  
}

#' spatial_binstat
#' 
#' creating binary matrix with specified feature (unavailable temporarily)
#'
#' @param spatial_adjust_obj list, 
#' @param bins vector, 
#' @param min.cutoff numeric, 
#' @param min.cnt numeric, 
#' @param pct.cutoff numeric, 
#' @param bin_border_col character, 
#' @param bin_target_col character, 
#' @param bin_border_lwd numeric, 
#' @param bin_target_lwd numeric, 
#' @param bin_border_lty numeric, 
#' @param bin_target_lty numeric, 
#' @param do.stat boolean, 
#' @param mtx_msk matrix, 
#' @param plot boolean, 
#' @param label boolean, 
#' @param pt_cex numeric
#'
#' @return list
#' 
spatial_binstat = function(spatial_adjust_obj, bins=c(8,8), min.cutoff=0, min.cnt=10, pct.cutoff=0.5,
                           bin_border_col='gray',bin_target_col='red',bin_border_lwd=0.5,bin_target_lwd=1,bin_border_lty=3,bin_target_lty=3,
                           do.stat=T, mtx_msk=NULL,
                           plot=T,label=F, pt_cex=0.5){
  
  new_img = spatial_adjust_obj$img_cord
  proj_name = spatial_adjust_obj$proj_name
  vec_col = spatial_adjust_obj$vec_col
  vec = spatial_adjust_obj$vec_feature
  feature = spatial_adjust_obj$feature
  
  if(is.null(min.cutoff)) min.cutoff = mean(vec)
  
  row_range = range(new_img[,2])
  row_indx = seq(row_range[1],row_range[2],length.out=bins[1]+1)
  col_range = range(new_img[,1])
  col_indx = seq(col_range[1],col_range[2],length.out=bins[2]+1)
  row_bin_center = (row_indx[2]-row_indx[1])/2
  col_bin_center = (col_indx[2]-col_indx[1])/2
  
  stat_mtx_pct = matrix(NA,nrow=bins[1],ncol=bins[2])
  stat_mtx_cnt = matrix(NA,nrow=bins[1],ncol=bins[2])
  stat_mtx_msk = matrix(0,nrow=bins[1],ncol=bins[2])
  if(plot){
    plot(0,type='n',ylim=range(new_img[,2]),xlim=range(new_img[,1]),main=feature, xaxt='n',yaxt='n',xlab='',ylab='')
    for(i in 1:nrow(new_img)){
      points(new_img[i,1],new_img[i,2],col=vec_col[i],pch=16,cex=pt_cex)
    }
    for(i in 1:(length(row_indx)-1)){
      for(j in 1:(length(col_indx)-1)){
        
        target_bins_indx = which(new_img[,2] >= row_indx[i] & new_img[,2] < row_indx[i+1] & new_img[,1] >= col_indx[j] & new_img[,1] < col_indx[j+1])
        if(length(target_bins_indx)>0){
          target_vec = vec[which(names(vec) %in% rownames(new_img)[target_bins_indx])]
          stat_mtx_pct[bins[1]-i+1,j] = ifelse(length(target_bins_indx) > min.cnt,length(which(target_vec>min.cutoff))/length(target_bins_indx),NA)
          stat_mtx_cnt[bins[1]-i+1,j] = length(target_bins_indx)
        }
        
        if(is.null(mtx_msk)){
          if(!is.na(stat_mtx_pct[bins[1]-i+1,j]) & stat_mtx_pct[bins[1]-i+1,j] > pct.cutoff){
            stat_mtx_msk[bins[1]-i+1,j] = 1
            rect(xleft = col_indx[j], ybottom = row_indx[i], xright = col_indx[j+1], ytop = row_indx[i+1], 
                 border = bin_target_col, lwd = bin_target_lwd, lty = bin_target_lty)
          }else{
            rect(xleft = col_indx[j], ybottom = row_indx[i], xright = col_indx[j+1], ytop = row_indx[i+1], 
                 border = bin_border_col, lwd = bin_border_lwd, lty = bin_border_lty)
          }
        }else{
          if(!is.na(mtx_msk[bins[1]-i+1,j]) & mtx_msk[bins[1]-i+1,j]==1){
            rect(xleft = col_indx[j], ybottom = row_indx[i], xright = col_indx[j+1], ytop = row_indx[i+1], 
                 border = bin_target_col, lwd = bin_target_lwd, lty = bin_target_lty)
          }else{
            rect(xleft = col_indx[j], ybottom = row_indx[i], xright = col_indx[j+1], ytop = row_indx[i+1], 
                 border = bin_border_col, lwd = bin_border_lwd, lty = bin_border_lty)
          }
        }
        
        if(label) text(col_indx[j]+col_bin_center,row_indx[i]+row_bin_center,round(stat_mtx_pct[bins[1]-i+1,j],2))
        
      }
    }
  }
  
  stat_summary = NULL
  stat_df = NULL
  if(do.stat){
    if(!is.null(mtx_msk)) stat_mtx_msk = mtx_msk
    stat_df = data.frame(feature=as.vector(stat_mtx_pct),group=as.vector(stat_mtx_msk))
    if(length(table(stat_df$group))==2 & min(table(stat_df$group))>1){
      stat_summary = list(
        'group_1' = mean(stat_df$feature[which(stat_df$group==1)],na.rm=T),
        'group_0' = mean(stat_df$feature[which(stat_df$group==0)],na.rm=T),
        'statistics' = wilcox.test(stat_df$feature[which(stat_df$group==1)],stat_df$feature[which(stat_df$group==0)])
      )
    }else{
      cat('[WARNING] too few groups for statistics:',table(stat_df$group),'..\n')
    }
  }
  
  spatial_adjust_obj$binstat = list(
    stat_mtx_pct = stat_mtx_pct,
    stat_mtx_cnt = stat_mtx_cnt,
    stat_mtx_msk = stat_mtx_msk,
    stat_summary = stat_summary,
    stat_df = stat_df,
    row_indx = row_indx,
    col_indx = col_indx,
    min.cutoff = min.cutoff,
    min.cnt = min.cnt,
    pct.cutoff = pct.cutoff
  )
  
  return(
    spatial_adjust_obj
  )
  
}

#' spatial_cordstat
#' 
#' analyzing Jaccard indexes between two spatial_adjust objects with specified features and visualizing spatial map 
#'
#' @param spatial_adjust_obj1 list, a spatial_adjust object
#' @param spatial_adjust_obj2 list, a spatial_adjust object with another specified feature name
#' @param min.cutoffs vector, a cutoff value on the feature to create a binary representation. If not specified (NULL), a mean value will be generated automatically. 
#' @param min.pct_cutoffs vector, a percentile cutoff value on the feature to create a binary representation. This parameter is mutually exclusive with min.cutoffs. If both are provided, min.pct_cutoffs takes precedence. Default, NULL. 
#' @param operator_steps vector, a vector containing step sizes for simulating the displacement of spatial object in eight cardinal directions. Default, c(4,4,4,4). 
#' @param do.stat boolean, determines whether the Jaccard indexes should be calculated. Default, TRUE. 
#' @param plot boolean, determines whether visulization should be generated. Default, TRUE. 
#' @param plot_bin boolean, determines whether visulization should be generated in binary style. Default, FALSE. 
#' @param plot_bin_col vector, a sequence of colors to use. Default, c('orange','blue'). 
#' @param cord_mov character, a specified displacement of spatial map for visualization. Valid values are 'add_X', 'minus_X', 'add_Y', 'minus_Y', 'add_X_add_Y', 'add_X_minus_Y', 'minus_X_add_Y', and 'minus_X_minus_Y'. Default, NULL. 
#' @param pt_cex numeric, size of point. Default, 0.5. 
#'
#' @return list
#' @export
spatial_cordstat = function(spatial_adjust_obj1,spatial_adjust_obj2, min.cutoffs=c(0,0),min.pct_cutoffs=NULL,
                            operator_steps=c(4,4,4,4), do.stat=T, plot=T, plot_bin=F, plot_bin_col=c('orange','blue'), cord_mov = NULL, pt_cex=0.5){
  
  new_img1 = spatial_adjust_obj1$img_cord
  proj_name1 = spatial_adjust_obj1$proj_name
  vec_col1 = spatial_adjust_obj1$vec_col
  vec1 = spatial_adjust_obj1$vec_feature
  feature1 = spatial_adjust_obj1$feature
  
  new_img2 = spatial_adjust_obj2$img_cord
  proj_name2 = spatial_adjust_obj2$proj_name
  vec_col2 = spatial_adjust_obj2$vec_col
  vec2 = spatial_adjust_obj2$vec_feature
  feature2 = spatial_adjust_obj2$feature
  
  if(is.null(min.cutoffs) & is.null(min.pct_cutoffs)) min.cutoffs = c(mean(vec1),mean(vec2))
  
  if(!is.null(min.pct_cutoffs)){
    min.cutoffs = c(quantile(vec1,min.pct_cutoffs[1]),quantile(vec2,min.pct_cutoffs[2]))
  }
  
  vec1_bin = ifelse(vec1>min.cutoffs[1],1,-1)
  vec2_bin = ifelse(vec2>min.cutoffs[2],1,-1)
  
  new_img2_mov = list()
  new_img2_mov[['add_X']] = t(t(new_img2) + c(operator_steps[1],0))
  new_img2_mov[['minus_X']] = t(t(new_img2) + c(-operator_steps[2],0))
  new_img2_mov[['add_Y']] = t(t(new_img2) + c(0,operator_steps[3]))
  new_img2_mov[['minus_Y']] = t(t(new_img2) + c(0,-operator_steps[4]))
  new_img2_mov[['add_X_add_Y']] = t(t(new_img2) + c(operator_steps[1],operator_steps[3]))
  new_img2_mov[['add_X_minus_Y']] = t(t(new_img2) + c(operator_steps[1],-operator_steps[4]))
  new_img2_mov[['minus_X_add_Y']] = t(t(new_img2) + c(-operator_steps[2],operator_steps[3]))
  new_img2_mov[['minus_X_minus_Y']] = t(t(new_img2) + c(-operator_steps[2],-operator_steps[4]))
  
  rownames(new_img1) = paste0('S',new_img1[,1],'_',new_img1[,2])
  rownames(new_img2) = paste0('S',new_img2[,1],'_',new_img2[,2])
  names(vec1_bin) = rownames(new_img1)
  names(vec2_bin) = rownames(new_img2)
  for(i in 1:length(new_img2_mov)){
    tmp = new_img2_mov[[i]]
    rownames(new_img2_mov[[i]]) = paste0('S',unlist(tmp[,1]),'_',unlist(tmp[,2]))
  }
  
  if(do.stat){
    new_img2_mov_stat = as.data.frame(t(sapply(1:length(new_img2_mov),function(i,y){
      x = new_img2_mov[[i]]
      names(y) = rownames(x)
      cm = intersect(rownames(x),rownames(new_img1))
      rs = c(
        round(length(which(vec1_bin[cm] + y[cm]>1))/(length(which(vec1_bin[cm]>0))+length(which(y[cm]>0))-length(which(vec1_bin[cm] + y[cm]>1))),4),
        round(length(which(vec1_bin[cm] + y[cm]>1))/max(c(length(which(vec1_bin[cm]>0)),length(which(y[cm]>0)))),4),
        round(length(which(vec1_bin[cm] + y[cm]>1))/min(c(length(which(vec1_bin[cm]>0)),length(which(y[cm]>0)))),4),
        round(length(which(vec1_bin + y>1))/(length(which(vec1_bin>0))+length(which(y>0))-length(which(vec1_bin + y>1))),4),
        round(length(which(vec1_bin + y>1))/max(c(length(which(vec1_bin>0)),length(which(y>0)))),4),
        round(length(which(vec1_bin + y>1))/min(c(length(which(vec1_bin>0)),length(which(y>0)))),4),
        length(which(vec1_bin[cm] + y[cm]>1)),
        length(which(vec1_bin>0)),
        length(which(y>0)),
        length(cm),
        length(vec1_bin)
      )
      names(rs) = c('pct_p','pct_pmin','pct_pmax','pct_op','pct_opmin','pct_opmax','cnt_p','cnt_p1','cnt_p2','cnt_px','cnt_all')
      rs
    },y=vec2_bin)))
    rownames(new_img2_mov_stat) = names(new_img2_mov)
    new_img2_mov_stat$diff_pct = new_img2_mov_stat$pct_p - new_img2_mov_stat$pct_op
  }
  
  if(plot){

    par(mar = c(5, 4, 4, 6) + 0.1, xpd = TRUE)
    
    row_range = range(new_img1[,2])
    col_range = range(new_img1[,1])
    new_img2_mov_tmp = new_img2
    if(!is.null(cord_mov)){
      new_img2_mov_tmp = new_img2_mov[[cord_mov]]
      row_range = range(c(row_range,range(new_img2_mov_tmp[,2])))
      col_range = range(c(col_range,range(new_img2_mov_tmp[,1])))
    }
    plot(0,type='n',ylim=row_range,xlim=col_range,main=paste0(feature1,':',feature2), xaxt='n',yaxt='n',xlab='',ylab='')
    
    if(plot_bin){
      vec_col1 = ifelse(vec1_bin>0,plot_bin_col[1],'white')
      vec_col2 = ifelse(vec2_bin>0,plot_bin_col[2],'white')
    }
    for(i in 1:nrow(new_img1)){ 
      points(new_img1[i,1],new_img1[i,2],col=vec_col1[i],pch=16,cex=pt_cex)
    }
    for(i in 1:nrow(new_img2_mov_tmp)){
      points(new_img2_mov_tmp[i,1],new_img2_mov_tmp[i,2],col=vec_col2[i],cex=pt_cex)
    }

    usr <- par("usr")
    legend(x = usr[2] + (usr[2]-usr[1])*0.05, bty='n',
       y = usr[4] - (usr[4]-usr[3])/3,   
       legend = c(feature1, feature2),
       col = c(vec_col1[1],vec_col2[1]),
       pch = c(16, 1))
    par(mar = c(5, 4, 4, 2) + 0.1, xpd = FALSE)
    
  }
  
  return(list(
    mov_summary = apply(new_img2_mov_stat,2,mean),
    mov_stat = new_img2_mov_stat,
    feature = paste0(feature1,':',feature2)
  ))
  
}
