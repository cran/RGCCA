## ---- setup, include=FALSE----------------------------------------------------
options(prompt = 'R> ', continue = '+ ')
options(ggrepel.max.overlaps = Inf)

def.chunk.hook  <- knitr::knit_hooks$get("chunk")
knitr::knit_hooks$set(chunk = function(x, options) {
  x <- def.chunk.hook(x, options)
  paste0("\n \\", "footnotesize","\n\n", x, "\n\n \\normalsize")
})

knitr::opts_chunk$set(
  fig.path = "figures/"
)

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("RGCCA")

## -----------------------------------------------------------------------------
RGCCA::available_methods()

## -----------------------------------------------------------------------------
library(RGCCA)
data(Russett)
colnames(Russett)

## -----------------------------------------------------------------------------
A <- list(
  Agric = Russett[, c("gini", "farm", "rent")],
  Ind = Russett[, c("gnpr", "labo")],
  Polit = Russett[, c("inst", "ecks",  "death", "demostab", "dictator")])

lab <- factor(
  apply(Russett[, 9:11], 1, which.max),
  labels = c("demost", "demoinst", "dict")
)

## -----------------------------------------------------------------------------
#Define the design matrix C.
C <- matrix(c(0, 0, 1,
              0, 0, 1,
              1, 1, 0), 3, 3)

C

## -----------------------------------------------------------------------------
fit <- rgcca(blocks = A, connection = C,
             tau = 1, ncomp = 2,
             scheme = "factorial",
             scale = TRUE,
             scale_block = FALSE,
             comp_orth = TRUE,
             verbose = FALSE)

## -----------------------------------------------------------------------------
print(fit)

## ---- fig.height = 12, fig.width=18, fig.cap = 'Block-weight vectors of a fitted RGCCA model.', fig.pos = "H"----
plot(fit, type = "weight", block = 1:3, comp = 1,
     display_order = FALSE, cex = 2)

## ---- fig.align='center', fig.cap = '\\label{fig:sample}Graphical display of the countries by drawing the block component of the first block against the block component of the second block, colored according to their political regime.', fig.height = 12, fig.width=18, fig.pos = "H"----
plot(fit, type = "sample",
     block = 1:2, comp = 1,
     resp = lab, repel = TRUE, cex = 2)

## ---- fig.align='center', fig.cap = 'Average Variance Explained of the different blocks.', fig.height = 8, fig.width=18, fig.pos = "H"----
plot(fit, type = "ave", cex = 2)

## ---- fig.align='center', fig.cap = 'Correlation circle associated with the first two components of the Agriculture block.', fig.height = 12, fig.width=18, fig.pos = "H"----
plot(fit, type = "cor_circle", block = 1, comp = 1:2, 
     display_blocks = 1:3, cex = 2)

## ---- fig.align='center', fig.cap = 'Biplot associated with the first two components ofthe Agriculture block.', fig.height = 12, fig.width=18, fig.pos = "H"----
plot(fit, type = "biplot", block = 1, 
     comp = 1:2, repel = TRUE, 
     resp = lab, cex = 2,
     show_arrow = TRUE)

## ---- cache = TRUE, message = FALSE-------------------------------------------
boot_out <- rgcca_bootstrap(fit, n_boot = 500, n_cores = 1)

## ---- size = "tiny"-----------------------------------------------------------
print(boot_out, block = 1:3, ncomp = 1)

## ---- fig.cap = 'Bootstrap confidence intervals for the block-weight vectors.', fig.height = 12, fig.width=18, fig.pos = "H"----
plot(boot_out, type = "weight", 
     block = 1:3, comp = 1, 
     display_order = FALSE, cex = 2,
     show_stars = TRUE)

## -----------------------------------------------------------------------------
fit.mcoa <- rgcca(blocks = A, method = "mcoa", ncomp = 2)

## -----------------------------------------------------------------------------
print(fit.mcoa)

## ---- fig.align='center', fig.cap = 'Biplot of the countries obtained by crossing the two first components of the superblock. Individuals are colored according to their political regime and variables according to their block membership.', fig.height = 12, fig.width=18, fig.pos = "H"----
plot(fit.mcoa, type = "biplot", 
     block = 4, comp = 1:2, 
     response = lab, 
     repel = TRUE, cex = 2)

## -----------------------------------------------------------------------------
fit <- rgcca(blocks = A, connection = C,
             tau = "optimal", scheme = "factorial")

## -----------------------------------------------------------------------------
fit$call$tau

## -----------------------------------------------------------------------------
set.seed(123)
perm_out <- rgcca_permutation(blocks = A, connection = C,
                              par_type = "tau",
                              par_length = 10,
                              n_cores = 1,
                              n_perms = 10)

## ---- width=30----------------------------------------------------------------
print(perm_out)

## ---- fig.height = 12, fig.width=18, fig.pos = "H", fig.cap = "Values of the objective function of RGCCA against the sets of tuning parameters, triangles correspond to evaluations on non-permuted datasets."----
plot(perm_out, cex = 2)

## -----------------------------------------------------------------------------
fit <- rgcca(perm_out)

## ----eval = FALSE-------------------------------------------------------------
#  # Download the dataset's package at http://biodev.cea.fr/sgcca/.
#  # --> gliomaData_0.4.tar.gz
#  if (!("gliomaData" %in% rownames(installed.packages()))) {
#    destfile <- tempfile()
#    download.file("http://biodev.cea.fr/sgcca/gliomaData_0.4.tar.gz", destfile)
#    install.packages(destfile, repos = NULL, type = "source")
#  }

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(eval = "gliomaData" %in% rownames(installed.packages()))

## -----------------------------------------------------------------------------
#  data(ge_cgh_locIGR, package = "gliomaData")
#  
#  blocks <- ge_cgh_locIGR$multiblocks
#  Loc <- factor(ge_cgh_locIGR$y)
#  levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#  blocks[[3]] <- Loc
#  
#  # check dimensions of the blocks
#  vapply(blocks, NCOL, FUN.VALUE = 1L)

## -----------------------------------------------------------------------------
#  fit.rgcca <- rgcca(blocks = blocks, response = 3, ncomp = 2, verbose = FALSE)

## -----------------------------------------------------------------------------
#  fit.rgcca$call$connection
#  fit.rgcca$call$tau

## -----------------------------------------------------------------------------
#  fit.rgcca$primal_dual

## -----------------------------------------------------------------------------
#  system.time(
#    rgcca(blocks = blocks, response = 3)
#  )

## ---- fig.align='center', fig.height = 12, fig.width=18, fig.pos = "H", fig.cap = "Graphical display of the tumors obtained by crossing the block components, and colored according to their location."----
#  plot(fit.rgcca, type = "sample", block = 1:2,
#       comp = 1, response = Loc, cex = 2)

## -----------------------------------------------------------------------------
#  fit.sgcca <- rgcca(blocks = blocks, response = 3, ncomp = 2,
#                     sparsity = c(0.0710, 0.2000, 1),
#                     verbose = FALSE)

## -----------------------------------------------------------------------------
#  print(fit.sgcca)

## ---- cache = TRUE------------------------------------------------------------
#  set.seed(27) #my favorite number
#  inTraining <- caret::createDataPartition(
#    blocks[[3]], p = .75, list = FALSE
#  )
#  training <- lapply(blocks, function(x) as.matrix(x)[inTraining, , drop = FALSE])
#  testing <- lapply(blocks, function(x) as.matrix(x)[-inTraining, , drop = FALSE])
#  
#  cv_out <- rgcca_cv(blocks = training, response = 3,
#                     par_type = "sparsity",
#                     par_value = c(.2, .2, 0),
#                     par_length = 10,
#                     prediction_model = "lda",
#                     validation = "kfold",
#                     k = 3, n_run = 5, metric = "Accuracy",
#                     n_cores = 2)

## -----------------------------------------------------------------------------
#  print(cv_out)

## ---- fig.height = 12, fig.width=18, fig.pos = "H", fig.cap = "Accuracies of the models on the different validation folds against the sets of tuning parameters."----
#  plot(cv_out, cex = 2)

## -----------------------------------------------------------------------------
#  fit <- rgcca(cv_out)
#  print(fit)

## -----------------------------------------------------------------------------
#  pred <- rgcca_predict(fit, blocks_test = testing, prediction_model = "lda")

## -----------------------------------------------------------------------------
#  pred$confusion$test

## -----------------------------------------------------------------------------
#  projection <- rgcca_transform(fit, blocks_test = testing)

## ---- cache = TRUE, message = FALSE-------------------------------------------
#  fit_stab <- rgcca_stability(fit,
#                              keep = vapply(
#                                fit$a, function(x) mean(x != 0),
#                                FUN.VALUE = 1.0
#                              ),
#                              n_boot = 100, verbose = TRUE, n_cores = 2)

## ---- fig.height = 12, fig.width=18, fig.align='center', fig.cap = 'Graphical display of the tumors obtained by crossing the components of GE1 and CGH1, and colored according to their location.', fig.pos = "H"----
#  plot(fit_stab, type = "sample", block = 1:2,
#       comp = 1, resp = as.character(Loc)[inTraining],
#       cex = 2
#       )

## ---- cache = TRUE------------------------------------------------------------
#  boot_out <- rgcca_bootstrap(fit_stab, n_boot = 500)

## ---- fig.height = 12, fig.width=18, fig.pos = "H", fig.cap = "Bootstrap confidence intervals for the block-weight vectors associated to block GE."----
#  plot(boot_out, block = 1,
#       display_order = FALSE,
#       n_mark = 50, cex = 1.5, cex_sub = 17,
#       show_star = TRUE)

