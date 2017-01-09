require(ape)
require(stringr)

parseTree <- function (control_file) {
  lines <- readLines(control_file)
  nwk   <- grep("^ +\\(.+$", lines, value=T)
  tree  <- read.tree(text=nwk)
  return(tree)
}

parseTheta <- function (control_file) {
  lines <- readLines(control_file)
  line  <- grep("thetaprior", lines, value=T)
  theta_prior <- regexec("thetaprior\\s+=\\s+(\\d+)\\s+(\\d+)", line)
  theta_prior <- as.numeric(regmatches(line, theta_prior)[[1]][-1])
  names(theta_prior) <- c("alpha", "beta")
  return(theta_prior)
}

parseTau <- function (control_file) {
  lines <- readLines(control_file)
  line  <- grep("tauprior", lines, value=T)
  tau_prior <- regexec("tauprior\\s+=\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)", line)
  tau_prior <- as.numeric(regmatches(line, tau_prior)[[1]][-1])
  names(tau_prior) <- c("alpha", "beta", "dirichlet")
  return(tau_prior)
}

getSettings <- function (wdir=".") {
  wdir <- gsub("/$", "", wdir)

  # How many repeats
  dirs <- list.dirs(wdir, recursive=F, full.names=F)
  repeats <- grep("^repeat\\d+$", dirs, value=T)
  
  # sets
  analyses  <- list.dirs(paste(wdir, repeats[1], sep="/"), recursive=F, full.names=F)
  splitList <- strsplit(analyses, "\\.")
  
  sets <- unique(unlist(lapply(splitList, "[", 5)))
  
  # priors
  prior_pairs <- unique(unlist(lapply(splitList, "[", 2)))
  prior_pairs_splits <- strsplit(prior_pairs, "-")
  names_thetas <- unique(unlist(lapply(prior_pairs_splits, "[", 1)))
  names_taus   <- unique(unlist(lapply(prior_pairs_splits, "[", 2)))
  
  thetas <- vector("list", length(names_thetas)); names(thetas) <- names_thetas
  taus   <- vector("list", length(names_taus)); names(taus) <- names_taus
  
  for (i in 1:length(names_thetas)) {
    if (sets[1] %in% c("mol", "prior", "mol-NNI")) {
      prog <- "BPP"
    } else {
      prog <- "iBPP"
    }
    analysis <- paste0("std.", names_thetas[i], "-tau1.Tree1.", prog, ".", sets[1])
    thetas[[i]] <- parseTheta(paste0(wdir, "/repeat0/", analysis,"/", analysis, ".ctl"))
  }
  
  for (i in 1:length(names_taus)) {
    if (sets[1] %in% c("mol", "prior", "mol-NNI")) {
      prog <- "BPP"
    } else {
      prog <- "iBPP"
    }
    analysis <- paste0("std.theta1-", names_taus[i], ".Tree1.", prog, ".", sets[1])
    taus[[i]] <- parseTau(paste0(wdir, "/repeat0/", analysis,"/", analysis, ".ctl"))
  }
  
  # trees
  TreeNames <- unique(unlist(lapply(splitList, "[", 3)))
  noTrees   <- length(TreeNames)
  
  trees <- vector("list", noTrees)
  for (i in 1:noTrees) {
    if (sets[1] %in% c("mol", "prior", "mol-NNI")) {
      prog <- "BPP"
    } else {
      prog <- "iBPP"
    }
    analysis <- paste("std", prior_pairs[1], TreeNames[i], prog, sets[1], sep=".")
    trees[[i]] <- parseTree(paste0(wdir, "/repeat0/", analysis,"/", analysis, ".ctl"))
  }
  names(trees) <- TreeNames
  return(list(wdir=wdir,
              analyses=analyses,
              noRepeats=length(repeats),
              repeats=repeats,
              trees=trees,
              prior_pairs=prior_pairs,
              tetha_priors=thetas,
              tau_priors=taus,
              sets=sets))
}


plotGuideTrees <- function (filename, settings, ...) {
  n <- length(settings$trees)
  pdf(filename, width=4*n, height=4, useDingbats=F)
  layout(matrix(1:n, nrow=1))
  par(mar=c(0,0,1,0), omi=c(0,0,0,0)+.1)
  for (i in 1:n) {
    main <- paste("Tree", i)
    tree <- settings$trees[[i]]
    tree <- makeNodeLabel(tree,prefix="split")
    plot.phylo(tree, show.node.label=T, main=main, ...)
  }
  dev.off()
}


# parse node support and ESS values from a BPP output file
parseBPP <- function (BPPoutfile, ESSwarning=T) {
  lines  <- readLines(BPPoutfile, warn=F)
  myline <- grep("Guide tree with posterior probability for presence of nodes", lines)
  if (length(myline)<1) stop("No output-tree found in ", BPPoutfile) # A meaningfull error message when this fails. (Sometimes there is no output in the .out.txt-file...)
  nwk    <- lines[myline+1]
  tree   <- read.tree(text=nwk)
  
  # get support values
  supports <- tree$node.label                      # get node values in ape-order
  if (is.null(supports)) {                         # this happens when BPP was stuck in the one species model
    supports <- rep(0, tree$Nnode)
  } else {
    supports[which(supports == "")] <- 0           # assign '0' to nodes without support (BPP omits the support value when it's 0)
    matches  <- gregexpr("(^'?#?)|('$)", supports) # clean them
    regmatches(supports, matches, invert=F) <- ""  # the variable 'supports' get cleaned here
  }
  names(supports) <- makeNodeLabel(tree, prefix="split")$node.label
  
  # get MAP tree
  MAPline <- grep("Summarizing the posterior of parameters under the MAP tree ", lines)
  MAP <- sub("Summarizing the posterior of parameters under the MAP tree ","", lines[MAPline])
  if (length(MAP) == 0)
    MAP = NA
  
  # get ESS values
  myline <- grep("ESS\\*", lines)
  if (length(myline) > 0) {
    myESS <- lines[myline]
    myESS <- strsplit(myESS, " +")[[1]][-1]
    myESS <- as.numeric(myESS)
  } else {
    myESS <- NA
    if (ESSwarning) {
      warning("ESS values not found", immediate.=T)
    }
  }
  return(list(supports=supports, MAP=MAP, ESS=myESS))
}


parseBPPunguided <- function (BPPoutfile, ESSwarning=F) {
  lines <- readLines(BPPoutfile, warn=F)
  lines <- gsub("(^\\s+)|(\\s+$)", "", lines)
  myLine_B_i <- grep("^\\(B\\)", lines)
  myLine_C_i <- grep("^\\(C\\)", lines)
  myLine_D_i <- grep("^\\(D\\)", lines)
  
  if (length(myLine_B_i)<1 || length(myLine_C_i)<1 || length(myLine_D_i)<1) {
    warning("BPP outfile incomplete: ", BPPoutfile)
    return(list(models=data.frame(models=NA, pp=NA), species=data.frame(sp=NA, pp=NA)))
  }

  # parse species delimitations
  sp_delims <- data.frame(model=vector("character"), pp=vector("numeric"), stringsAsFactors=F)
  i <- 1
  for (line_i in (myLine_B_i+2):(myLine_C_i-2)) {
    sp_delims[i,] <- str_match(lines[line_i], "^\\d+\\s+([0-9.]+)\\s+\\d+\\s+\\((.+)\\)$")[,3:2]
    i <- i+1
  }
  
  # parse species pp
  sp_pp <- data.frame(sp=vector("character"), pp=vector("numeric"), stringsAsFactors=F)
  i <- 1
  for (line_i in (myLine_C_i+2):(myLine_D_i-2)) {
    sp_pp[i,] <- str_match(lines[line_i], "^\\d+\\s+([0-9.]+)\\s+(.+)$")[,3:2]
    i <- i+1
  }
  
  if (ESSwarning) {
    warning("Currently there are no ESS values output by BPP...", immediate.=T)
  }
  
  return(list(models=sp_delims, species=sp_pp))
}


list_means <- function (list) {
  n <- length(list[[1]]) # number of repeats
  v <- vector("numeric")
  for (i in 1:n) {
    values <- lapply(list, "[", i)
    values <- unlist(lapply(list, "[", i))
    v[i] <- median(values)
  }
  vv <- matrix(v, nrow=nrow(list[[1]]))
  rownames(vv) <- rownames(list[[1]])
  colnames(vv) <- colnames(list[[1]])
  return(vv)
}


get_all_scenarios <- function (list) {  # list is a raw list of unguided analyses results (UGout in parseBPPmulti)
  nrep <- length(list)      # number of repeats
  nanl <- length(list[[1]]) # number of analyses
  
  # first get all models and all species that were sampled
  models  <- vector("character")
  species <- vector("character")
  for (i in 1:nrep) {
    for (j in 1:nanl) {
      models <- c(models, list[[i]][[j]]$models$model)
      species <- c(species, list[[i]][[j]]$species$sp)
    }
  }
  
  models  <- unique(na.omit(models))
  species <- unique(na.omit(species))
  
  return(list(models=models, species=species))
}


means_unguided <- function (list, settings, models, species, limit=0.05) { # list is a raw list of unguided analyses results (UGout in parseBPPmulti)

  # now get the pps for all models and species
  # These must become lists with trees as elements
  # column names of the tables must be prior combinations (so better loop through analyses here and best also above...)
  models_pp  <- list()
  species_pp <- list()
  models_total_means  <- data.frame(row.names=models)
  species_total_means <- data.frame(row.names=species)
  for (tree in names(settings$trees)) {
    models_pp[[tree]]  <- data.frame(row.names=models)
    species_pp[[tree]] <- data.frame(row.names=species)
    for (pair in settings$prior_pairs) {
      analysis <- paste("std", pair, tree, "BPP.mol-NNI", sep=".")
      models_temp  <- matrix(NA, nrow=length(models),  ncol=settings$noRepeats, dimnames=list(models=models))
      species_temp <- matrix(NA, nrow=length(species), ncol=settings$noRepeats, dimnames=list(models=species))
      for (i in 1:settings$noRepeats) {
        pp <- list[[i]][[analysis]]$models$pp; names(pp) <- list[[i]][[analysis]]$models$model
        models_temp[,i] <- as.numeric(pp[models])
        pp <- list[[i]][[analysis]]$species$pp; names(pp) <- list[[i]][[analysis]]$species$sp
        species_temp[,i] <- as.numeric(pp[species])
      }
      # replace all NA with 0 (those modes or species were never sampled)
      models_temp[which(is.na(models_temp), arr.ind=T)] <- 0
      species_temp[which(is.na(species_temp), arr.ind=T)] <- 0
      
      # get median of all repeats
      models_pp[[tree]][[pair]]  <- models_temp  <- apply(models_temp, 1, median)
      species_pp[[tree]][[pair]] <- species_temp <- apply(species_temp, 1, median)
    }
    # total means over all pairs (for sorting)
    models_total_means[[tree]]  <- apply(models_pp[[tree]],  1, mean)
    species_total_means[[tree]] <- apply(species_pp[[tree]], 1, mean)
  }
  
  # sorting and filtering
  models_total_means[["total"]]  <- apply(models_total_means,  1, mean)
  species_total_means[["total"]] <- apply(species_total_means, 1, mean)

  models_total_means  <- models_total_means[order(models_total_means$total, decreasing=T),]
  species_total_means <- species_total_means[order(species_total_means$total, decreasing=T),]
  
  models_total_means  <- models_total_means[which(models_total_means$total >= limit),]
  species_total_means <- species_total_means[which(species_total_means$total >= limit),]
  
  for (tree in names(settings$trees)) {
    models_pp[[tree]]  <- models_pp[[tree]][rownames(models_total_means),]
    species_pp[[tree]] <- species_pp[[tree]][rownames(species_total_means),]
  }
  
  return(list(models=models_pp, species=species_pp))
}


parseBPPmulti <- function (settings, unguided_limit=0.05) {
  odir <- getwd()
  setwd(settings$wdir)
  
  # output will be a list of large tables:
  # the items of the list are the repeats
  # the tables have split numbers as columns and the analysis as rownames
  
  out <- list()
  MAPtrees <- list()
  UGout <- list()
  ESS <- list()
  
  for (rep in settings$repeats) {
    cat(rep,"...\n")
    # parse GUIDED analyses
    guided_analyses <- grep("NNI", settings$analyses, value=T, invert=T)
    nr <- length(guided_analyses)
    nc <- settings$trees$Tree1$Nnode
    nodelabel <- paste("split", 1:nc, sep="")
    table <- matrix(nrow=nr, ncol=nc, dimnames=list(guided_analyses, nodelabel))
    ESSvals <- list()
    MAPvals <- matrix(NA, nrow=length(settings$analyses), ncol=nc, dimnames=list(settings$analyses, paste0("split",1:nc)))
    for (analysis in guided_analyses) {
      BPPoutfile <- paste(rep, "/", analysis, "/", analysis, ".out.txt", sep="")
      BPPoutput <- parseBPP(BPPoutfile, ESSwarning=F)
      table[analysis,] <- as.numeric(BPPoutput$supports[nodelabel])
      ESSvals[[analysis]] <- BPPoutput$ESS
      if (is.na(BPPoutput$MAP)) {
        MAPvals[analysis,] <- rep(NA, nc)
      } else {
        MAPvals[analysis,] <- as.numeric(strsplit(BPPoutput$MAP, "")[[1]])
      }
    }
    out[[rep]] <- table
    ESS[[rep]] <- ESSvals
    MAPtrees[[rep]] <- MAPvals
    
    # parse UNGUIDED analyses
    unguided_analyses <- grep("NNI", settings$analyses, value=T, invert=F)
    for (analysis in unguided_analyses) {
      BPPoutfile <- paste(rep, "/", analysis, "/", analysis, ".out.txt", sep="")
      UGout[[rep]][[analysis]] <- parseBPPunguided(BPPoutfile)
    }
  }
  
  # medians of guided analyses
  mean_table <- list_means(out)
  
  
  # get all species models and all species that were sampled in the unguided analyses
  if (length(UGout) == 0) {
    UG_tables <- NA
    UGmean_tables <- NA
  } else {
    scenarios <- get_all_scenarios(UGout)
    models  <- scenarios$models
    species <- scenarios$species

  
    # medians of unguided analyses
    UGmean_tables <- means_unguided(UGout, settings, models, species, limit=unguided_limit)
    
    # make readable tables of single repeats
    UG_tables <- list()
    for (rep in settings$repeats) {
      UG_tables[[rep]] <- list()
      for (tree in names(settings$trees)) {
        UG_tables[[rep]][["models"]] <- list()
        UG_tables[[rep]][["species"]] <- list()
        UG_tables[[rep]][["models"]][[tree]]  <- data.frame(row.names=models)
        UG_tables[[rep]][["species"]][[tree]] <- data.frame(row.names=species)
        for (pair in settings$prior_pairs) {
          analysis <- paste("std", pair, tree, "BPP.mol-NNI", sep=".")
          
          pp <- as.numeric(UGout[[rep]][[analysis]]$models$pp)
          names(pp) <- UGout[[rep]][[analysis]]$models$model
          UG_tables[[rep]][["models"]][[tree]][[pair]] <- pp[models]
          UG_tables[[rep]][["models"]][[tree]][[pair]][which(is.na(UG_tables[[rep]][["models"]][[tree]][[pair]]))] <- 0
          
          pp <- as.numeric(UGout[[rep]][[analysis]]$species$pp)
          names(pp) <- UGout[[rep]][[analysis]]$species$sp
          UG_tables[[rep]][["species"]][[tree]][[pair]] <- pp[species]
          UG_tables[[rep]][["species"]][[tree]][[pair]][which(is.na(UG_tables[[rep]][["species"]][[tree]][[pair]]))] <- 0
        }
        # AT THIS POINT IT WOULD BE EASY TO DETECT ANALYSES THAT WERE STUCK IN THE ONE SPECIES MODEL
        # sort and filter like UG_mean-tables
        UG_tables[[rep]][["models"]][[tree]]  <- UG_tables[[rep]][["models"]][[tree]][rownames(UGmean_tables$models[[tree]]),]
        UG_tables[[rep]][["species"]][[tree]] <- UG_tables[[rep]][["species"]][[tree]][rownames(UGmean_tables$species[[tree]]),]
      }
    }
  }
  
  setwd(odir)
  
  return(list(supports=out,
              mean_supports=mean_table,
              MAPtrees=MAPtrees,
              unguided=UG_tables,
              unguided_mean_supports=UGmean_tables,
              ESS=ESS))
}


extractSet <- function (BPPresultsMatrix, tree, set) {
  if (set=="mol-NNI") {
    out <- list()
    out[["models"]]  <- BPPresultsMatrix$models[[tree]]
    out[["species"]] <- BPPresultsMatrix$species[[tree]]
  } else {
    pattern <- paste(tree, ".+", set, sep="")
    ii <- grep(pattern, rownames(BPPresultsMatrix))
    out <- BPPresultsMatrix[ii,]
  }
  return(out)
}


write.BPP.tables <- function (dir, settings, BPPresults) {
  dir <- gsub("(^\\s*.?/)|(/\\s*$)", "", dir)
  if (!dir.exists(dir)) dir.create(dir)
  
  for (tree in names(settings$trees)) {
    for (set in settings$sets) {
      if (set == "mol-NNI") {
        table_models  <- extractSet(BPPresults$unguided_mean_supports, tree, set)$models
        table_species <- extractSet(BPPresults$unguided_mean_supports, tree, set)$species
        file <- paste(dir, "/", tree, "_", set, "_means.tsv", sep="")
        write.table(table_models, file, append=F, quote=F, sep="\t", row.names=T, col.names=NA)
        cat("Posterior probabilities of delimited species\n", file=file, append=T)
        write.table(table_species, file, append=T, quote=F, sep="\t", row.names=T, col.names=F)
        for (rep in settings$repeats) {
          table_models  <- extractSet(BPPresults$unguided[[rep]], tree, set)$models
          table_species <- extractSet(BPPresults$unguided[[rep]], tree, set)$species
          file <- paste(dir, "/", tree, "_", set, rep, ".tsv", sep="")
          write.table(table_models, file, append=F, quote=F, sep="\t", row.names=T, col.names=NA)
          cat("Posterior probabilities of delimited species\n", file=file, append=T)
          write.table(table_species, file, append=T, quote=F, sep="\t", row.names=T, col.names=F)
        }
      } else {
        table <- extractSet(BPPresults$mean_supports, tree, set)
        write.table(table, paste(dir, "/", tree, "_", set, "_means.tsv", sep=""), append=F, quote=F, sep="\t", row.names=T, col.names=NA)
        for (rep in settings$repeats) {
          table <- extractSet(BPPresults$supports[[rep]], tree, set)
          write.table(table, paste(dir, "/", tree, "_", set, "_", rep, ".tsv", sep=""), append=F, quote=F, sep="\t", row.names=T, col.names=NA)
        }
      }
    }
  }
}


get.colors <- function (x, colors, stops) {
  my.colors <- c()
  for (number in x) {
    i=1
    numbers.color <- "black"
    for (col in colors) {
      if (stops[i] <= number & number <= stops[i+1]) {
        numbers.color <- col
        break
      } else {
        i=i+1
      }
    }
    my.colors <- c(my.colors, numbers.color)
  }
  return(my.colors)
}


plot.pp <- function (tree, table, nrow=3, ncol=3, colors=c("tomato", "orange1", "lightgoldenrod1", "lightgreen"), colstops=c(0, .3, .6, .9, 1), sep=c(.2,.2), exp=1, add.order=F, ...) {
  n <- nrow(table)
  tree <- makeNodeLabel(tree, prefix="split")
  
  adjx <- rep(scale(seq(from=0, by=sep[1], length.out=ncol), center=T, scale=F)+.5, nrow)
  adjy <- rev(rep(scale(seq(from=0, by=sep[2], length.out=nrow), center=T, scale=F)+.5, each=ncol))
  
  adjxt <- rev(rep(scale(seq(from=0, by=sep[1]*10, length.out=ncol), center=T, scale=F)+.5, nrow))
  adjyt <- rep(scale(seq(from=0, by=sep[2]*10, length.out=nrow), center=T, scale=F)+.5, each=ncol)
  
  plot.phylo(tree, label.offset=.2, no.margin=T, ...)
  i = 1
  for (pair in rownames(table)) {
    x <- table[pair,]
    x <- x[tree$node.label] # This is probably not necessary
    col <- get.colors(x, colors, colstops)
    nodelabels(pch=22, col="gray30", lwd=.3, bg=col, cex=exp, adj=c(adjx[i], adjy[i]))
    if (add.order) {
      nodelabels(i, frame="none", adj=c(adjxt[i], adjyt[i]), cex=exp/4.25)
      cat(i, ": ", pair, "\n", sep="")
    }
    i = i+1
  }
}


plot.all.pp <- function (settings, BPPresults, order=c("prior", "mol", "trait", "intgr"), colors=c("tomato", "orange1", "lightgoldenrod1", "lightgreen"), colstops=c(0, .3, .6, .9, 1), ...) {
  sets <- settings$sets
  sets <- grep("NNI", sets, value=T, invert=T) # exclude unguided analyses
  sets <- sort(factor(sets, levels=order, ordered=T))
  cat("Order of plots from left to right:", names(settings$trees), "\n")
  cat("Order of plots from top to bottom:", as.character(sets), "\n")
  nTrees <- length(settings$trees)
  nSets  <- length(sets)
  layout(matrix(1:(nTrees*nSets), nrow=nSets))
  for (tree in names(settings$trees)) {
    phylo <- settings$trees[[tree]]
    for (set in sets) {
      table <- extractSet(BPPresults$mean_supports, tree, set)
      plot.pp(phylo, table, colors=colors, colstops=colstops, ...)
    }
  }
  layout(1)
}


# write.publication.tables <- function (dir=getwd(), tables) {
#   file <- paste(dir, "/publication_tables.tsv", sep="")
#   cat("Table S12. Posterior probabilities of iBPP analyses, nicely sorted.\n\n\n", file=file, append=F)
#   types <- c("BPP.prior", "BPP.mol", "iBPP.trait", "iBPP.intgr")
#   prior.combinations <- c("theta1.tau1", "theta2.tau1", "theta3.tau1", "theta1.tau2", "theta2.tau2", "theta3.tau2", "theta1.tau3", "theta2.tau3", "theta3.tau3")
#   repeats <- paste("repeat", 0:9, sep="")
#   for (tree in names(trees)) {  
#     for (type in types) {
#       for (prior.combination in prior.combinations) {
#         analysis <- paste("std", prior.combination, tree, type, sep=".")
#         table <- data.frame(row.names=paste("split", 1:7, sep="")) # create table to save pps of all repeats and splits
#         for (rep in repeats) {
#           table[[rep]] <- tables[[rep]][[analysis]][1:7] # also discard convergence row
#         }
#         # sort splits
#         table <- table[trees[[tree]]$node.label,]
#         cat(tree, " - ", type, " - ", prior.combination, "\n", file=file, sep="", append=T)
#         write.table(table, file, append=T, quote=F, sep="\t", row.names=T, col.names=NA)
#         cat("\n\n", file=file, sep="", append=T)
#       }
#     }
#   }
# }
