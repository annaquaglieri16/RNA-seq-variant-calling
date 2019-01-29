buildKEGGIdx <- function (entrezIDs, species = "human", min.size = 1, updated = FALSE, 
                          exclude = c()) 
{
  species = normalizeSpecies(species)
  updatedSuccess = FALSE
  entrezIDs = as.character(entrezIDs)
  kegg = NULL
  print("Building KEGG pathways annotation object ... ")
  if (!updated) {
    kegg.pathways = loadKeggData()
    kegg = kegg.pathways[[species.fullToShort[[tolower(species)]]]]
  }
  else {
    kegg = kegg.gsets(species = species.fullToShort[[tolower(species)]], 
                 id.type = "kegg")
      updatedSuccess = TRUE
  }
  if (is.null(kegg)) 
    stop("Failed to load the KEGG pathway data.")
  gsets = kegg$kg.sets
  gsets.ez = gsets
  gsets = ids2indices(gsets.ez, entrezIDs, remove.empty = TRUE)
  gsets.ez = gsets.ez[names(gsets)]
  gsets.ids = sapply(names(gsets.ez), function(x) as.character(substr(x, 
                                                                      1, 8)))
  gsets.names = sapply(names(gsets.ez), function(x) as.character(substr(x, 
                                                                        10, nchar(x))))
  tmp = rep(NA, length(kegg$kg.sets))
  tmp[kegg$sig.idx] = "Signaling"
  tmp[kegg$met.idx] = "Metabolism"
  tmp[kegg$dise.idx] = "Disease"
  anno = data.frame(ID = gsets.ids, GeneSet = gsets.names, 
                    NumGenes = paste0(sapply(gsets, length), "/", sapply(gsets.ez, 
                                                                         length)), Type = tmp[match(names(gsets), names(kegg$kg.sets))])
  rownames(anno) = gsets.names
  names(gsets) = gsets.names
  names(gsets.ez) = gsets.names
  if (!updatedSuccess) {
    db.info = egsea.data(simple = TRUE, returnInfo = TRUE)
    ver = db.info$kegg$info$version
    dat = db.info$kegg$info$date
  }
  else {
    ver = "NA"
    dat = date()
  }
  gs.annot = GSCollectionIndex(original = gsets.ez, idx = gsets, 
                               anno = anno, featureIDs = entrezIDs, species = species, 
                               name = "KEGG Pathways", label = "kegg", version = ver, 
                               date = dat)
  gs.annot = selectGeneSets(gs.annot, min.size = min.size)
  if (length(gs.annot@idx) == 0) 
    cat("KEGG pathway collection is empty.\\n")
  if (length(exclude) > 0) {
    sel = !tolower(gs.annot@anno[, "Type"]) %in% tolower(exclude)
    gs.annot@idx = gs.annot@idx[sel]
    gs.annot@original = gs.annot@original[sel]
    gs.annot@anno = gs.annot@anno[sel, ]
  }
  return(gs.annot)
}


normalizeSpecies <- function(species){
  human.names = c("human", "homo sapiens", "hs")
  mouse.names = c("mouse", "mus musculus", "mm")
  rat.names = c("rat", "rattus norvegicus" , "rn")
  species = tolower(species)
  if (species %in% human.names){
    species = "Homo sapiens"
  }else if (species %in% mouse.names){
    species = "Mus musculus"
  }
  else if (species %in% rat.names){
    species = "Rattus norvegicus"
  }
  else{
    stop("Unrecognized species.")
  }
  return(species)
}

species.fullToShort = list()
species.fullToShort[["homo sapiens"]] = "human"
species.fullToShort[["mus musculus"]] = "mouse"
species.fullToShort[["rattus norvegicus"]] = "rat"

GSCollectionIndex <- setClass(
  "GSCollectionIndex",            
  slots = c(original = "list",
            idx = "list",
            anno = "data.frame",                
            featureIDs = "character",
            species = "character",
            name = "character",
            label = "character",
            version="character",
            date="character"),                      
  prototype = list(original = list(),
                   idx = list(),
                   anno = data.frame(),                
                   featureIDs = "",
                   species = "",
                   name = "",
                   label = "",
                   version="",
                   date = "")           
)

setGeneric(name="selectGeneSets",
           def = function(object, gs.names=NULL, min.size=1){
             standardGeneric("selectGeneSets")
           }
)

setMethod(f = "selectGeneSets",
          signature(object = "GSCollectionIndex"),
          definition = function(object, gs.names=NULL, min.size=1){
            if (length(object@idx) == 0)
              return(object)
            gs.annot.top = GSCollectionIndex()
            if (is.null(gs.names)){
              gs.names = names(object@idx[sapply(object@idx, function(x) 
                length(x)) >= min.size])    
            }   
            else{
              gs.names = gs.names[sapply(object@idx[gs.names], function(x) 
                length(x)) >= min.size]
            }
            if ("GeneSet" %in% colnames(object@anno)){
              sel = match(gs.names, object@anno[, "GeneSet"])
              gs.annot.top@original = object@original[sel]
              gs.annot.top@idx = object@idx[sel]
              gs.annot.top@anno = object@anno[sel,]
              gs.annot.top@anno = droplevels(gs.annot.top@anno)
            }else{
              warning("The 'GeneSet' column was not found in 'anno'.")
              gs.annot.top@original = object@original
              gs.annot.top@idx = object@idx
            }
            gs.annot.top@label = object@label
            gs.annot.top@featureIDs = object@featureIDs
            gs.annot.top@species = object@species
            gs.annot.top@name = object@name
            gs.annot.top@version = object@version
            gs.annot.top@date = object@date
            return(gs.annot.top)          
          }
          
)


  