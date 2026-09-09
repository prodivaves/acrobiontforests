Ecological niche comparison scripts
This folder contains R scripts associated with the ecological niche analyses in:

Sánchez-González, L. A., Prieto-Torres, D. A., & Espinosa-Chávez, O. J. (2026). Ecological Niche Differences Underlie the Assembly of Endemic Birds in Acrobiont Forests of Northern Mesoamerica. Journal of Biogeography, 53, e70141. https://doi.org/10.1111/jbi.70141

Scripts
01_select_variables_spearman_vif.R
Implements the first environmental-variable selection approach described in the article:

extract values from the 19 WorldClim 1.4 bioclimatic variables at occurrence localities;
calculate a Spearman correlation matrix;
remove strongly correlated predictors using usdm::vifcor(th = 0.8, method = "spearman");
remove remaining multicollinear predictors using usdm::vifstep(th = 10);
save the retained variables, correlation matrices, VIF tables, and session information.
The article also evaluated a second PCA-derived six-variable approach. That procedure is not reconstructed in this script because the exact original PCA-selection rule is not contained in the supplied VIF script. If the original PCA-selection code is available, it should be archived as a separate script rather than inferred retrospectively.

02_niche_overlap_equivalency_similarity.R
Performs pairwise PCA-env niche comparisons using ecospat:

common environmental-variable screening using Spearman |r| < 0.8 followed by VIF < 10;
PCA-env transformation;
Schoener's D niche overlap;
niche equivalency test;
directional background similarity tests;
1,000 randomizations for equivalency and similarity tests.
The humboldt package is intentionally omitted from this repository version, as requested.

All paths in the scripts are relative to the repository root. Alternatively, define the environment variable NICHE_PROJECT_ROOT with the absolute path to the repository before running the scripts.
