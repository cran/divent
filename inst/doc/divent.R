## ----global_options, include=FALSE--------------------------------------------
set.seed(97310)

## ----load_paracou6------------------------------------------------------------
library("divent")
paracou_6_abd
# Number of individuals in each community
abd_sum(paracou_6_abd)

## ----plot_paracou6------------------------------------------------------------
autoplot(paracou_6_abd[1, ])

## ----rcommunity---------------------------------------------------------------
rc <- rcommunity(1, size = 10000, distribution = "lnorm")
autoplot(rc, fit_rac = TRUE, distribution = "lnorm")

## ----estimation---------------------------------------------------------------
div_richness(paracou_6_abd)
ent_shannon(paracou_6_abd)
ent_simpson(paracou_6_abd)

## ----naive_shannon------------------------------------------------------------
library("dplyr")
paracou_6_abd %>% 
  as_probabilities() %>% 
  ent_shannon()

## ----shannon_estimators-------------------------------------------------------
ent_shannon(paracou_6_abd)
ent_shannon(paracou_6_abd, estimator = "ChaoJost")

## ----ent_tsallis--------------------------------------------------------------
ent_tsallis(paracou_6_abd, q = 1)

## ----div_hill-----------------------------------------------------------------
div_hill(paracou_6_abd, q = 1)

## ----lnq----------------------------------------------------------------------
(d2 <- div_hill(paracou_6_abd, q = 2)$diversity)
ln_q(d2, q = 2)
(e2 <- ent_tsallis(paracou_6_abd, q = 2)$entropy)
exp_q(e2, q = 2)

## ----PhyloDiversity-----------------------------------------------------------
div_phylo(paracou_6_abd, tree = paracou_6_taxo, q = 1)

## ----DivVector----------------------------------------------------------------
# Richness of a community of 100 species, each of them with 10 individuals
div_richness(rep(10, 100))

## ----SBDiversity--------------------------------------------------------------
# Similarity is computed from the functional distance matrix of Paracou species
Z <- fun_similarity(paracou_6_fundist)
# Calculate diversity of order 2
div_similarity(paracou_6_abd, similarities = Z, q = 2)

## ----div_profile--------------------------------------------------------------
profile_hill(paracou_6_abd) %>% autoplot

## ----PDiversityProfile--------------------------------------------------------
profile_phylo(paracou_6_abd, tree = paracou_6_taxo) %>% autoplot
# Similarity matrix
Z <- fun_similarity(paracou_6_fundist)
profile_similarity(paracou_6_abd, similarities = Z) %>% autoplot

## ----div_level----------------------------------------------------------------
# Estimate the diversity of 1000 individuals
div_hill(paracou_6_abd, q = 1, level = 1000)

## ----div_coverage-------------------------------------------------------------
# Estimate the diversity at 80% coverage
div_hill(paracou_6_abd, q = 1, level = 0.8)

## ----accum_hill---------------------------------------------------------------
accum_hill(
  paracou_6_abd[1, ], 
  q = 1, 
  levels = 1:500,
  n_simulations = 100
) %>% 
  autoplot()

## ----accum_div_phylo----------------------------------------------------------
accum_div_phylo(
  paracou_6_abd[1, ],
  tree = paracou_6_taxo,
  q = 1, 
  levels = 1:2000
) %>% 
  autoplot()

## ----MetaCommunitydf----------------------------------------------------------
# Abundances of three communities with four species
(abd <- matrix(
  c(
    10,  0, 25, 10, 
    20, 15, 10, 35, 
     0, 10,  5,  2
  ),
  ncol = 4
))
# Community weights
w <- c(1, 2, 1)

## -----------------------------------------------------------------------------
(communities <- as_abundances(abd, weights = w))

## ----MetaCommunityMC----------------------------------------------------------
(mc <- metacommunity(communities))
plot(communities, type = "Metacommunity")

## ----DivPart------------------------------------------------------------------
div_part(paracou_6_abd, q = 1)

## ----gamma--------------------------------------------------------------------
div_hill(paracou_6_abd, q = 1, gamma = TRUE)

