#' @include Algorithm.SDM.R
#' @import methods
#' @import SpaDES.tools
#' @importFrom raster raster stack extract predict reclassify layerStats calc
NULL

setGeneric("project", function(obj, Env, ...) {
  return(standardGeneric("project"))
})

setMethod("project", "Algorithm.SDM", function(obj, Env, ...) {
  model = get_model(obj, ...)
  ##Dividir el Modelo
  ###Crear la carpeta temporal
  #path <- get("tmpdir", envir = .PkgEnv)
  path <- tempdir()
  if (!(file.path(path, ".rasters_2") %in% list.dirs(path)))
  {  dir.create(paste0(path, "/.rasters_2"))}
  
  n <- 2
  
  for (i in seq_len(length(Env@layers))) {
     splitRaster(Env[[i]], n, n, path = paste0(path, "/.rasters_2/", names(Env[[i]])))
  }
  ##
  proj_2 <- list()
  for(i in seq_len(n*n)){
  ##reclasificar
    listado_completo <- list.files(paste0(path, "/.rasters_2") , pattern = paste("*tile",i,".gri$", sep = ''), recursive = TRUE, full.names = TRUE )
    Env_temp <- raster::stack(listado_completo)
    proj_3 = suppressWarnings(raster::predict(Env_temp, model, fun = function(model,x) {
      x = as.data.frame(x)
      for (j in seq_len(length(Env_temp@layers))) {
        if (Env_temp[[j]]@data@isfactor) {
          x[, j] = as.factor(x[, j])
          x[, j] = droplevels(x[, j])
          levels(x[, j]) = Env_temp[[j]]@data@attributes[[1]]$ID
        }
      }
      return(predict(model, x))
    }))
    # Rescaling projection
    proj_2[i] = reclassify(proj_3, c(-Inf, 0, 0))
  }  
  ##Mosaic -> recomponiendo los rasters
  proj <- mergeRaster(proj_2)
  unlink(file.path(path, ".rasters_2"), recursive=TRUE)
  ##
  # Rescaling projection
  #proj = reclassify(proj, c(-Inf, 0, 0))
  if(all(obj@data$Presence %in% c(0,1))) # MEMs should not be rescaled
    proj = proj / proj@data@max
  names(proj) = "Projection"
  obj@projection = proj
  if(all(obj@data$Presence %in% c(0,1))) # MEMs can't produce binary
  obj@binary <- reclassify(proj, c(-Inf,obj@evaluation$threshold,0,
                                   obj@evaluation$threshold,Inf,1))
  gc()
  return(obj)

})

setMethod("project", "MAXENT.SDM", function(obj, Env, ...) {
  model = get_model(obj, Env, ...)
  ##Para devidir los rasters  de ENV
  path <- tempdir()
  if (!(file.path(path, ".rasters_2") %in% list.dirs(path)))
  {  dir.create(paste0(path, "/.rasters_2"))}
  
  n <- 2
  
  for (i in seq_len(length(Env@layers))) {
    splitRaster(Env[[i]], n, n, path = paste0(path, "/.rasters_2/", names(Env[[i]])))
  }
  ##
  proj_2 <- list()
  for(i in seq_len(n*n)){
    ##reclasificar
    listado_completo <- list.files(paste0(path, "/.rasters_2") , pattern = paste("*tile",i,".gri$", sep = ''), recursive = TRUE, full.names = TRUE )
    Env_temp <- raster::stack(listado_completo)
    
    proj_3 = raster::predict(Env_temp, model, fun = function(model, x) {
      x = as.data.frame(x)
      for (j in seq_len(length(Env@layers))) {
        if (Env[[j]]@data@isfactor) {
          x[, j] = as.factor(x[, j])
          x[, j] = droplevels(x[, j])
          levels(x[, j]) = Env_temp[[j]]@data@attributes[[1]]$ID
        }
      }
      return(predict(model, x))
    })
    # Rescaling projection
    proj_2[i] = reclassify(proj_3, c(-Inf, 0, 0))
  }
  ##Mosaic -> recomponiendo los rasters
  proj <- mergeRaster(proj_2)
  unlink(file.path(path, ".rasters_2"), recursive=TRUE)
  #proj = reclassify(proj, c(-Inf, 0, 0))
  ##
  if(!all(obj@data$Presence %in% c(0,1))) # MEMs should not be rescaled
    proj = proj / proj@data@max
  names(proj) = "Projection"
  obj@projection = proj
  if(all(obj@data$Presence %in% c(0,1))) # MEMs can't produce binary
    obj@binary <- reclassify(proj, c(-Inf,obj@evaluation$threshold,0,
                                     obj@evaluation$threshold,Inf,1))
  gc()
  return(obj)
})
