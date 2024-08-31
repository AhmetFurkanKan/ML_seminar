# ML_seminar
Seminar in Machine Learning in Econometrics
# Economic Predictions with Big Data: An Extension

This repository contains code used for the seminar paper, based on the original work by Domenico Giannone, Michele Lenza, and Giorgio E. Primiceri in their paper *"Economic Predictions With Big Data: The Illusion of Sparsity"*.

## Acknowledgements

The original code and methodology were developed by Domenico Giannone, Michele Lenza, and Giorgio E. Primiceri. Their work can be found [here](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA17842). This repository includes additional analysis and modifications made to their code for educational and research purposes.

## Modifications

The following changes were made to the original code:
•	SpikeSlabGLP_Lap8.m and SpikeSlabGLP_t2.m - The main function with Laplace and student's t likelihood for the posterior evaluation of the model, which can be applied to any dataset with a response variable y and a set of predictors x.
•	macro1_Laplace.m  and macro1_studentst.m - Matlab codes to evaluate the posterior distribution of the unknown parameters of the model in macro1 application. Furthermore it creates joint density distributions of the hyperparameters q and log(γ), and posterior density of q.
•	data - Datasets used 
