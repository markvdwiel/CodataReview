This repository supports the manuscript: 'Guiding adaptive shrinkage by co-data to improve regression-based prediction and feature selection' by Mark van de Wiel and Wessel van Wieringen.
Repository contains the following R-scripts: 

1. GroupLassoVSGroupAdaptiveLasso.R: scripts that reproduces simulation-based comparison between sparse group lasso and group-adaptive lasso
2. Xtune_Varbvs_merge.R: scripts that illustrates how to merge the functionalities of R-packages xtune and varbvs to estimate co-data guided feature-specific prior inclusion probabilities for the spike-and-slab model
3. squeezy_hyperDemo: demo script on how to perform hyperparameter shrinkage for the group-adaptive lasso
4. squeezy_hyper.R: source file for the Demo.
5. Demo_ecpc_mirdata.R: data example on miRNA data with ecpc

Moreover, it contains the data and co-data for the data example: XYcodata.Rdata
