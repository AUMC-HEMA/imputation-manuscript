merge_mod <- function(FCS1, FCS2, BBmarkers=NULL, th = length(BBmarkers)*0.20, normalize=TRUE, leftout=FALSE) {
  
  name1 = FCS1@name
  name2 = FCS2@name
  
  FCS1 <- FCS1[,sort(FCS1@markers)]
  FCS2 <- FCS2[,sort(FCS2@markers)]
  
  cat("=== CytoBackBone ===\n")
  
  if(normalize==TRUE){
    cat(paste0("Normalizing profiles\n"))
    normalizedProfiles <- normalizeQuantile(list(FCS1,FCS2),BBmarkers)
    FCS1               <- normalizedProfiles[[1]]
    FCS2               <- normalizedProfiles[[2]]
  }
  
  FCS1@markers[FCS1@markers %in% BBmarkers]        <- paste0("BB_",FCS1@markers[FCS1@markers %in% BBmarkers])
  FCS2@markers[FCS2@markers %in% BBmarkers]        <- paste0("BB_",FCS2@markers[FCS2@markers %in% BBmarkers])
  FCS1markersnotbb                                 <- FCS1@markers[!(FCS1@markers %in% paste0("BB_",BBmarkers))]
  FCS2markersnotbb                                 <- FCS2@markers[!(FCS2@markers %in% paste0("BB_",BBmarkers))]
  
  if(length(intersect(FCS1markersnotbb,FCS2markersnotbb))>0){
    message("----------ERROR----------")
    message("backbone markers are not consistent between the two cytometric profiles")
    message("backbone markers: ",appendLF=FALSE)
    message(paste0(BBmarkers,collapse=","))
    message("specific markers for FCS1: ",appendLF=FALSE)
    message(paste0(FCS1markersnotbb,collapse=","))
    message("specific markers for FCS2: ",appendLF=FALSE)
    message(paste0(FCS2markersnotbb,collapse=","))
    message("intersection: ",appendLF=FALSE)
    intersect = intersect(FCS1markersnotbb,FCS2markersnotbb)
    message(paste0(intersect,collapse=","))
    message("----------ERROR----------")
    stop()
  }
  
  FCS1markersnotbb                                 <- FCS1markersnotbb[FCS1markersnotbb %in% FCS2markersnotbb]
  FCS2markersnotbb                                 <- FCS2markersnotbb[FCS2markersnotbb %in% FCS1markersnotbb]
  FCS1@markers[FCS1@markers %in% FCS1markersnotbb] <- paste0(FCS1@markers[FCS1@markers %in% FCS1markersnotbb],".1")
  FCS2@markers[FCS2@markers %in% FCS2markersnotbb] <- paste0(FCS2@markers[FCS2@markers %in% FCS2markersnotbb],".2")
  
  FCS1                        <- FCS1[, FCS1@markers[grep("cluster",FCS1@markers,invert=TRUE)]]
  FCS2                        <- FCS2[, FCS2@markers[grep("cluster",FCS2@markers,invert=TRUE)]]
  cat(paste0("profile 1 contains ",format(FCS1@cell.nb,big.mark=",")," cells\n"))
  cat(paste0("profile 2 contains ",format(FCS2@cell.nb,big.mark=",")," cells\n"))
  cat("===\n")
  FCS1_BB                     <- FCS1[,FCS1@markers[grep("BB_",FCS1@markers)]]@intensities
  colnames(FCS1_BB)           <- FCS1@markers[grep("BB_",FCS1@markers)]
  FCS2_BB                     <- FCS2[,FCS2@markers[grep("BB_",FCS2@markers)]]@intensities
  colnames(FCS2_BB)           <- FCS2@markers[grep("BB_",FCS2@markers)]
  dist                        <- FNN::knnx.dist(FCS2_BB,FCS1_BB,k=1,algorithm="kd_tree")
  FCS1_excluded               <- NULL
  
  # Keep track of excluded rows
  FCS1_excl_row               <- NULL
  FCS2_excl_row               <- NULL

  if(sum(dist>th)!=0){
    FCS1_excluded             <- FCS1[dist>th]
    FCS1_excl_row             <- dist>th
  }
  FCS1                        <- FCS1[dist<th]
  
  dist                        <- FNN::knnx.dist(FCS1_BB,FCS2_BB,k=1,algorithm="kd_tree")
  FCS2_excluded               <- NULL
  if(sum(dist>th)!=0){
    FCS2_excluded             <- FCS2[dist>th]
    FCS2_excl_row             <- dist>th
  }
  FCS2                        <- FCS2[dist<th]
  cat(paste0("profile 1 has ",format(FCS1@cell.nb,big.mark=",")," cells can be potentialy matched\n"))
  cat(paste0("profile 2 has ",format(FCS2@cell.nb,big.mark=",")," cells can be potentialy matched\n"))
  cat("===\n")
  
  max <- min(FCS1@cell.nb,FCS2@cell.nb)
  cat(paste0("maximum of number of cells that can be matched by CytoBackBone = ",format(max,big.mark=",")),"\n")
  cat("===\n")
  
  ##############################################################################
  all_iterations <- list()
  n <- 0
  ##############################################################################
  
  FCS1_BB                       <- FCS1[,FCS1@markers[grep("BB_",FCS1@markers)]]@intensities
  colnames(FCS1_BB)             <- FCS1@markers[grep("BB_",FCS1@markers)]
  FCS2_BB                       <- FCS2[,FCS2@markers[grep("BB_",FCS2@markers)]]@intensities
  colnames(FCS2_BB)             <- FCS2@markers[grep("BB_",FCS2@markers)]
  knnx                          <- FNN::get.knnx(FCS2_BB,FCS1_BB,k=1,algorithm="kd_tree")
  idx                           <- knnx$nn.index
  dist                          <- knnx$nn.dist
  table_FCS1_to_FCS2            <- cbind(1:nrow(idx),idx)
  table_FCS1_to_FCS2[dist>th,2] <- -1
  colnames(table_FCS1_to_FCS2)  <- c("FCS1","FCS2")
  knnx                          <- FNN::get.knnx(FCS1_BB,FCS2_BB,k=1,algorithm="kd_tree")
  idx                           <- knnx$nn.index
  dist                          <- knnx$nn.dist
  table_FCS2_to_FCS1            <- cbind(1:nrow(idx),idx)
  table_FCS2_to_FCS1[dist>th,2] <- -1
  colnames(table_FCS2_to_FCS1)  <- c("FCS2","FCS1")
  idx                           <- paste0(table_FCS1_to_FCS2[,1],"-",table_FCS1_to_FCS2[,2]) %in% paste0(table_FCS2_to_FCS1[,2],"-",table_FCS2_to_FCS1[,1])
  idx_a                         <- table_FCS1_to_FCS2[idx,"FCS1"]
  idx_b                         <- table_FCS1_to_FCS2[idx,"FCS2"]
  data                          <- cbind(FCS1@intensities[idx_a,],FCS2@intensities[idx_b,])
  colnames(data)                <- c(paste0(FCS1@markers,".FCS1"),paste0(FCS2@markers,".FCS2"))
  cat(paste0("step #",1,": ",format(nrow(data),big.mark=",")," cells matched (",format(max-nrow(data),big.mark=",")," cells remaining)\n"))
  
  ##############################################################################
  ids <- cbind(idx_a, idx_b)
  n <- n + 1
  iteration <- rep(n, nrow(data))
  info <- cbind(ids, iteration)
  combined <- cbind(data, info)
  all_iterations[[n]] <- data.frame(combined)
  ##############################################################################
  
  if(max-nrow(data)!=0){	
    i <- 2
    while(TRUE){
      
      FCS1                          <- FCS1[-idx_a]
      FCS2                          <- FCS2[-idx_b]
      FCS1_BB                       <- FCS1[,FCS1@markers[grep("BB_",FCS1@markers)]]@intensities
      colnames(FCS1_BB)             <- FCS1@markers[grep("BB_",FCS1@markers)]
      FCS2_BB                       <- FCS2[,FCS2@markers[grep("BB_",FCS2@markers)]]@intensities
      colnames(FCS2_BB)             <- FCS2@markers[grep("BB_",FCS2@markers)]
      knnx                          <- FNN::get.knnx(FCS2_BB,FCS1_BB,k=1,algorithm="kd_tree")
      idx                           <- knnx$nn.index
      dist                          <- knnx$nn.dist
      table_FCS1_to_FCS2            <- cbind(1:nrow(idx),idx)
      table_FCS1_to_FCS2[dist>th,2] <- -1
      colnames(table_FCS1_to_FCS2)  <- c("FCS1","FCS2")
      knnx                          <- FNN::get.knnx(FCS1_BB,FCS2_BB,k=1,algorithm="kd_tree")
      idx                           <- knnx$nn.index
      dist                          <- knnx$nn.dist
      table_FCS2_to_FCS1            <- cbind(1:nrow(idx),idx)
      table_FCS2_to_FCS1[dist>th,2] <- -1
      colnames(table_FCS2_to_FCS1)  <- c("FCS2","FCS1")
      idx                           <- paste0(table_FCS1_to_FCS2[,1],"-",table_FCS1_to_FCS2[,2]) %in% paste0(table_FCS2_to_FCS1[,2],"-",table_FCS2_to_FCS1[,1])
      idx_a                         <- table_FCS1_to_FCS2[idx,"FCS1"]
      idx_b                         <- table_FCS1_to_FCS2[idx,"FCS2"]
      
      if(sum(idx)==0)
        break
      if(sum(idx)==1)
        data_s                    <- cbind(t(FCS1@intensities[idx_a,]),t(FCS2@intensities[idx_b,]))
      else
        data_s                    <- cbind(FCS1@intensities[idx_a,],FCS2@intensities[idx_b,])
      
      data                          <- rbind(data,data_s)
      
      ##############################################################################
      ids <- cbind(idx_a, idx_b)
      n <- n + 1
      iteration <- rep(n, nrow(data_s))
      info <- cbind(ids, iteration)
      colnames(data_s) <- c(paste0(FCS1@markers,".FCS1"),paste0(FCS2@markers,".FCS2"))
      combined <- cbind(data_s, info)
      all_iterations[[n]] <- data.frame(combined)
      ##############################################################################
      
      cat(paste0("step #",i,": ",format(nrow(data),big.mark=",")," cells matched (",format(max-nrow(data),big.mark=",")," cells remaining)","\n"))
      i <- i+1
      
      if((max-nrow(data))==0)
        break
      
    }
  }
  
  if(!is.null(FCS1_excluded)){
    FCS1_excluded          <- c(FCS1_excluded,FCS1)
    FCS1_excluded@markers  <- gsub("BB_","",FCS1_excluded@markers)
    FCS1_excluded@name     <- paste0("cells specific to ",name1)
  }
  
  if(!is.null(FCS2_excluded)){
    FCS2_excluded          <- c(FCS2_excluded,FCS2)
    FCS2_excluded@markers  <- gsub("BB_","",FCS2_excluded@markers)
    FCS2_excluded@name     <- paste0("cells specific to ",name2)
  }
  
  excluded = NULL
  if(!is.null(FCS1_excluded) && !is.null(FCS2_excluded)){
    excluded               <- c(FCS1_excluded,FCS2_excluded)
    excluded@markers       <- gsub("BB_","",excluded@markers)
    excluded@name          <- paste0("cells specific to ",name1," or ",name2)
  }
  
  FCS                    <- as.FCS(data)
  FCS                    <- FCS[,sort(FCS@markers)]
  FCS                    <- merge_header(FCS)
  FCS@markers            <- gsub("BB_","",FCS@markers)
  FCS@name               <- paste0(name1," + ",name2)
  
  cat("====================\n")
  
  if(leftout==FALSE){
    return(list('FCS'=FCS, 'metadata'=t <- do.call('rbind', all_iterations), 
                'FCS1_excluded'=FCS1_excl_row, 'FCS2_excluded'=FCS2_excl_row))
  }else{
    return(list(merged=FCS,specific.FCS1=FCS1_excluded,specific.FCS2=FCS2_excluded,specific=excluded))
  }
}



#' @title FCS class definition
#'
#' @description FCS are S4 objects containing marker cell expressions obtained from cytometry profiles.
#'
#' @details This object mainly stores for each cell the intensities of all cell markers.
#'
#' The slot 'trans.para' is a named list contains different parameters depending of the transformation applied on the marker expression intensities. The scale (cofactor) of the arcsinh transformation function is parametrized using the 'arcsinh.scale' value. The shift of the log transformation function is parametrized using the 'log.shift' value and the base of the log transformation function is parametrized using the 'log.base' value. If no transformation function have been applied, the 'trans.para' slot is set to NULL.
#'
#' @slot name a character indicating the internal name of the FCS object
#' @slot profiles a character vector containing the names of the FCS profiles 
#' @slot cell.nb an integer value indicating the number of FCS profiles
#' @slot markers a character vector containing the marker names
#' @slot markers.nb an integer value indicating the number of markers
#' @slot intensities a numeric matrix containing the intensities of each marker for each FCS profile
#' @slot trans a character specifying the name of a transformation function applied on the marker expression intensities. Possible values are "arcsinh" for arc sin hyperbolic transformation, "log" for logarithmic transformation, or "none" for no transformation
#' @slot trans.para a named list containing parameters of the transformation. Please refer to the details section for more details
#' @slot trans.exclude a character vector containing the marker names for which no transformation has been applied on
#'
#' @import methods
#'
#' @name FCS-class
#' @rdname FCS-class
#' @exportClass FCS
FCS <- setClass("FCS",
                slots=c(name          = "character",
                        profiles          = "character",
                        cell.nb           = "integer",
                        markers           = "character",
                        markers.nb        = "integer",
                        intensities       = "matrix",
                        trans             = "character",
                        trans.para        = "ANY",
                        trans.exclude     = "ANY"),
                validity=function(object){
                  if(length(object@profiles)!=length(unique(object@profiles)))
                    stop("Error in profiles slot: profile names are not unique")
                  if(length(object@markers)!=length(unique(object@markers)))
                    stop("Error in markers slot: marker names are not unique")
                  if(object@cell.nb!=length(object@profiles))
                    stop("Error in cell.nb slot: cell.nb do not correspond to the number of profile names")
                  if(object@markers.nb!=length(object@markers))
                    stop("Error in markers.nb slot: markers.nb do not correspond to the number of marker names")
                  if(!is.element(object@trans,c("arcsinh","log","none")))
                    stop("Error in trans slot: trans do not contain allowed value (allowed value are \"arcsinh\", \"log\", \"none\")")
                  if(!is.null(object@trans.para) && class(object@trans.para)[1]!="list")
                    stop("Error in trans.para slot: trans.para must be a of type list or NULL")
                  if(!is.null(object@trans.exclude) && class(object@trans.exclude)[1]!="character")
                    stop("Error in trans.exclude slot: trans.exclude must be a of type character or NULL")
                  return(TRUE)
                }
)
setMethod("initialize",c("FCS"),
          function(.Object,
                   name              = "",
                   profiles          = "",
                   cell.nb           = 0,
                   markers           = "",
                   markers.nb        = 0,
                   intensities       = as.matrix(0)){
            if(cell.nb==0)
              stop("Error can not create a FCS object with no profile")
            .Object@name              = name     
            .Object@profiles          = profiles
            .Object@cell.nb           = cell.nb
            .Object@markers           = markers
            .Object@markers.nb        = markers.nb
            .Object@trans             = "none"
            .Object@trans.para        = NULL
            .Object@trans.exclude     = NULL
            .Object@intensities       = intensities
            methods::validObject(.Object)
            return(.Object)
          }
)



#' @title Extraction of subsets of data from FCS objects
#'
#' @description Extracts a subset of data from a FCS object.
#'
#' @details For The parameter i represents a vector of cells to extract and the parameter j represents a vector of markers to extract.
#'
#' @param x a FCS object
#' @param i a numeric, logical or character vector
#' @param j a numeric, logical or character vector
#'
#' @return a S4 object of class FCS
#'
#' @name extract
#' @rdname extract-methods
NULL

#' @rdname extract-methods
#' @export
setMethod("[",c("FCS","ANY","ANY"),
          function(x,i,j){
            
            if(!missing(i) && length(i)>x@cell.nb)
              stop("Too many cluster profiles to extract")
            if(!missing(j) && length(j)>length(x@markers))
              stop("Too many markers to extract")
            if(!missing(i) && !(is.logical(i) || is.numeric(i) || is.character(i)))
              stop(paste0("Wrong types in i: ",typeof(i)))
            if(!missing(j) && !(is.logical(j) || is.numeric(j) || is.character(j)))
              stop(paste0("Wrong types in j: ",typeof(j)))
            
            if(missing(i)){
              k <- 1:x@cell.nb
            }else{
              if(is.character(i)){
                k <- which(x@profiles %in% i)
              }else{
                k <- i
              }
            }
            
            if(missing(j)){
              l <- 1:length(x@markers)
            }else{
              if(is.character(j)){
                l <- unlist(sapply(paste0("^",j,"$"),grep,x=x@markers))
              }else{
                l <- j
              }    
            }
            
            FCS <- FCS(name        = x@name,
                       profiles           = x@profiles[k],
                       cell.nb            = length(x@profiles[k]),
                       markers            = x@markers[l],
                       markers.nb         = length(x@markers[l]),
                       intensities        = x@intensities[k,l,drop=FALSE])
            
            return(FCS)
          }
)


# Combination of FCS objects
#
# @description Combines two or several FCS objects.
#
# @param x a first FCS object
# @param ... further FCS objects to be combined
#
# @return a S4 object of class FCS
#
# @name c
# @rdname c-methods
NULL

# @rdname c-methods
# @export
setMethod("c",c("FCS"),
          function(x,...){
            
            other.FCS  <- list(x,...)
            name        <- c()
            markers     <- x@markers
            cell.nb <- 0
            
            i <- 1
            for(FCS in other.FCS){
              if(class(FCS) != "FCS")
                stop(paste0("Cannot combine objects of different classes (element at position ",i," is of type ",class(FCS)))
              name        <- c(name,FCS@name)
              markers     <- union(markers,FCS@markers)
              cell.nb <- cell.nb+FCS@cell.nb
              i <- i+1
            }
            name        <- paste0(name,collapse=";")
            profiles    <- as.character(1:cell.nb)
            
            intensities <- matrix(NA,ncol=length(markers),nrow=cell.nb,dimnames=list(profiles,markers))
            
            i <- 1
            for(FCS in other.FCS){
              nb                                   <- FCS@cell.nb
              intensities[i:(i+nb-1),FCS@markers] <- FCS@intensities
              i <- i+nb
            }
            dimnames(intensities) <- NULL
            
            FCS <- FCS(name = name,
                       profiles      = profiles,
                       cell.nb   = length(profiles),
                       markers       = markers,
                       markers.nb    = length(markers),
                       intensities   = intensities)
            return(FCS)
          }
)


# @title Coercion to a FCS object
#
# @description Coerces a numeric matrix into a FCS object.
#
# This function transforms a numeric matrix into one FCS object.
#
# @details The matrix must have its column names corresponding to the cell markers.
#
# @param object a numeric matrix
# @param name a character specifying the internal name of the FCS object to create
#
# @return a S4 object of class FCS
#
# @name as.FCS
# @rdname as.FCS-methods
#
# @export
setGeneric("as.FCS", function(object,name="FCS") { standardGeneric("as.FCS") })

# @rdname as.FCS-methods
# @export
setMethod("as.FCS",c("matrix"),
          function(object,name){  
            data           <- object
            dimnames(data) <- NULL
            FCS <- FCS(name = name,
                       profiles      = as.character(1:nrow(object)),
                       cell.nb       = nrow(object),
                       markers       = colnames(object),
                       markers.nb    = ncol(object),
                       intensities   = data)
            return(FCS)
          }
)


# title Internal - Extract a dictionary from a FCS file
#
# @description This function is used internally to extract the correspondence between the original marker names (first column) and the true marker names (second column)
#
# @param fcs.file a character indicating the location of the fcs file containing the correspondences
#
# @return a two-column data.frame providing the correspondence between the original marker names (first column) and the true marker names (second column)
extract.dictionary <- function(fcs.file){
  flowframe       <- flowCore::read.FCS(fcs.file[1],trans = FALSE)
  dictionary      <- flowframe@parameters@data[, c(1, 2)]
  dictionary[, 1] <- make.names(dictionary[, 1])
  return(dictionary)
}


# title Internal - Renaming cell markers
#
# @description This function is used internally to rename the cell markers based on a dictionary.
#
# @details dictionary is a data.frame used to rename the marker names. The first column must correspond to the original marker names, the second column must correspond to the new marker names. 
#
# @param header a character vector containing the original maker names
# @param dictionary a character vector containing a correspondence between the original and the new marker names
#
# @return a character vector containing the renamed marker names
rename.markers <- function(header,dictionary){
  header         <- make.names(header)
  dictionary[,1] <- as.vector(dictionary[,1])
  dictionary[,2] <- as.vector(dictionary[,2])
  if(length(unique(dictionary[,1]))!=length(dictionary[,1])){
    stop("Duplicate in dictionary 'original marker names'")
  }
  occurences <- table(dictionary[,2])
  occurences <- occurences[occurences>1]
  redondant.names <- names(occurences)
  for (i in 1:nrow(dictionary)) {
    if (any(dictionary[i, 2] %in% redondant.names)) {
      temp             <- gsub("X\\.(.*)\\.Di(_clust$|$)","(\\1)Di\\2",dictionary[i,1])
      dictionary[i,2] <- paste0(dictionary[i,2],"-",temp)
    }
    if(!is.na(dictionary[i,2])){
      header[which(header == dictionary[i,1])[1]] <- dictionary[i,2]
    }
  }
  return(header)
}


# title Internal - Removing of cell markers to exclude from a matrix
#
# @description This function is used internally to remove one or several cell markers from a numeric matrix.
#
# @param data a numeric matrix
# @param exclude a character vector containing the cell markers to be excluded
#
# @return a numeric matrix without the cell markers to exclude
exclude.markers <- function(data,exclude){
  exclude.flags <- exclude %in% colnames(data)
  if(any(!(exclude.flags))){
    warning(paste0("Unknown marker to exclude: ",paste(exclude[!exclude.flags],collapse=", ")))
  }
  data    <- data[,!(colnames(data) %in% exclude)]
  return(data)
}


merge_header <- function(FCS){
  bb              <- FCS@intensities[,grep("BB_",FCS@markers)]
  colnames(bb)    <- FCS@markers[grep("BB_",FCS@markers)]
  res <- c()
  M   <- c()
  for(i in seq(1,ncol(bb),by=2)){
    tmp <- apply(cbind(bb[,i],bb[,i+1]),1,mean)
    res <- cbind(res,tmp)
    M   <- c(M,colnames(bb)[i])
  }
  colnames(res) <- gsub(".FCS1","",M)
  evo             <- FCS@intensities[,grep("BB_",FCS@markers,invert=TRUE)]
  colnames(evo)   <- FCS@markers[grep("BB_",FCS@markers,invert=TRUE)]
  colnames(evo)   <- gsub(".FCS1","",colnames(evo))
  colnames(evo)   <- gsub(".FCS2","",colnames(evo))
  FCS  <- cbind(res,evo)
  FCS <- as.FCS(FCS)
  return(FCS)
}
