# Modeling and Predicting Spatio-temporal Dynamics of PM2.5 Concentrations Through Time-evolving Covariance Models

Welcome to the repository for the implementation of the paper "Modeling and Predicting Spatio-temporal Dynamics of PM2.5 Concentrations Through Time-evolving Covariance Models" by Ghulam A. Qadir and Ying Sun.

The paper is available on [arXiv](https://arxiv.org/abs/2202.12121).

## Overview

This repository contains the R code implementation of spatio-temporal dynamics modeling and prediction of PM2.5 concentrations using time-evolving covariance models. The code utilizes various R packages for statistical and spatial analysis.

## Requirements

To run the code in this repository, you need to have R installed along with the following packages:

- `MASS`
- `mvtnorm`
- `fields`
- `doParallel`
- `scoringRules`

You can install these packages in R using the following commands:

```r
install.packages(c("MASS", "mvtnorm", "fields", "doParallel", "scoringRules"))
