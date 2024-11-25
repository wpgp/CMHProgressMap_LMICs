# CMHProgressMap_LMICs
 
25-11-2024

------------------------------------------------------------------------

**Table 1.** Files and their descriptions within the CMHProgressMap_LMICs GitHub
repository 

| Name       | Type   | Description                                                                                                                                                                                                                                                                                                                                                                                                  |
|------------|--------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| covariates | Folder | This folder contains a demo of the format of the data extracted from geospatial covariates considered when modelling the health and development indicators. This file is required to run all R scripts in this repository. Examples of geospatial covariate datasets can be found from <https://hub.worldpop.org/project/categories?id=14>.                                                                     |
| indicators | Folder | This folder contains a demo of subnational maternal and child health indicators to model. This file is required to run all R scripts in this repository. The indicators were extracted from the 2014 (DHS-7)[1] and 2022 (DHS-8) [2] for Kenya , which are publicly available after registration onto the Measure DHS website (www.dhsprogram.com).                     |
| raster     | Folder | This folder contains the file structure for the geospatial covariates considered when modelling the health and development indicators. This file is required to run all R scripts in this repository. Examples of geospatial covariate datasets can be found from <https://hub.worldpop.org/project/categories?id=14>.                 |
| results    | Folder | Folder to contain all prediction and uncertainty gridded datasets (raster files), county and subcounty aggregations and model information produced from the workflow file                                                                                                                                                                                     |
| shapes     | Folder | This folder contains the shapefiles required to run all R scripts in this repository. These should be the administrative boundaries of the study area as polygons and the location of the clusters in the study area as points (lat/lon) These shapefiles can be obtained from the DHS program at www.dhsprogram.com.                                                                                        |
| raster_preparation_covariates.R | R    | This file contains code for the creation of the covariates from the raster files                     |
| workflow   | R      | R script for modelling the subnational maternal and child health indicators. The files required to run this script are the covariates and indicator csv files and the files in the shp folder.                                                        |
| utilities  | R      | R script containing helper functions for the main workflow file 

------------------------------------------------------------------------

# Script for modelling the child and maternal health indicators - workflow.R

This code was prepared for modelling indicators in Kenya however it has been presented 
in a general way to allow application to other contexts.

First covariates are then scaled and centred based on the associated rasters.
The geospatial covariate selection is two-staged.In the first stage, we check
for multicollinearity amongst the geospatial covariates.In the second stage,
we employ the backward stepwise model selection method.

To check for multicollinearity, a Pearson correlation matrix for the
geospatial covariates is created and any pairs with a Pearson
correlation coefficient r>0.8 are flagged. The flagged covariates are then individually fitted in
non-Bayesian binomial generalised linear models (GLMs). The Bayesian
information criteria (BIC) of the models are then calculated. The
covariate in the model with a lower BIC is retained while the covariate
in the model with the greater BIC is omitted for the target indicator.

After checking for multicollinearity, a backward model selection
algorithm is used to select the best (sub)set of geospatial covariates
for the target indicator. The algorithm is as follows. The remaining
geospatial covariates are fitted in a non-Bayesian binomial GLM and the
BIC is calculated. A covariate is removed from the model and the BIC is
recalculated. If the recalculated BIC is less than the previously
calculated BIC, this subset of covariates is preferred. These steps are
performed iteratively until the recalculated BIC is not less than the
BIC calculated from the previous iteration. At this point, the best
(sub)set of geospatial covariates have been attained and they will be
used when constructing the Bayesian point-referenced spatial binomial
GLM in INLA.

A Bayesian point-referenced spatial binomial GLM is then estimated, using R-INLA
the model contains a random effect with Matern covariance based on a mesh 
constructed from the country boundary and the GPS points where the indicators 
are collected, random effects corresponding to the district and state 
in which the point is located are also included [4,5]. Fixed effects based 
on the covariates identified in the model selection outlined are also included. 
More information on this model can be found in
https://www.nature.com/articles/s41597-023-01961-2 or for a wider introduction 
to INLA please see https://www.paulamoraga.com/book-geospatial/sec-inla.html. 

Additional components must be constructed before fitting the model in
INLA. First a mesh of the study domain is constructed with the shape
file and coordinates within the target indicator file. Using this mesh
object, a stochastic partial differential equation (SPDE) object is
defined with functions in INLA where the priors of the spatial decay
parameter and spatial variance parameter is defined. With the mesh
object, INLA stack “A” matrices are created and stacked with the INLA
stack functions. Finally, these components, along with the model are
fitted into the INLA function.

Once the model is fitted 1000 posterior samples are taken from the
INLA model. The script then reads in the
raster files corresponding to the geospatial covariates of the model for
the target indicators and compiles it as a prediction data frame.
Finally, the predicted values are computed from the prediction data
frame, INLA mesh objects, INLA posterior sample objects and the random 
effects corresponding to the district and state. This results in 1000 
samples from the modelled distribution for each 1x1km grid square.

This matrix of samples is produced for the indicator both round 1 corresponding to the 2014 DHS survey [1]
and round 2 corresponding to the 2022 DHS survey[2]. A matrix of 
samples for the change in the indicator is then produced for the change in the 
indicator between the two surveys. 

These matrices of samples at the 1x1km level are then aggregated to both 
subcounty and county administrative boundaries using population weighted aggregation
using population data is from Worldpop unconstrained population raster 1km 
resolution [6]. The mean, standard deviation 95% credible intervals 
are produced directly from the 1000 samples. In the case of the aggregation 
of the change in the indicator a further column is included corresponding to 
the probability that the true change is greater than 0. 
The aggregated ouputs are availible at https://doi.org/10.5258/SOTON/WP00791

The script then outputs raster tif files corresponding to the mean, median, standard 
deviation and 95% credible interval of the indicators at the 5x5km level. 
These high resolution outputs are availible at https://doi.org/10.5258/SOTON/WP00790




------------------------------------------------------------------------

# Acknowledgement

The authors acknowledge the support of the PMO Team at WorldPop.
Moreover, the authors would like to thank the DHS Program staff for their input on the construction of
some of the indicators. This work was approved by the ethics and research governance committee at the 
University of Southampton.

# Suggested citation
Dorey, P.,Chan, H.M.T., Tejedor-Garavito, N., Priyatikanto, R., Bonnie, A., Williams, E.M., Johnson, M., Tatem, A.J., Pezzulo, C. 2024. CMHProgressMap_LMICs: Child and Maternal Health Indicator Mapping from DHS Data in LMICs, version 1.1. WorldPop, University of Southampton. 

# References

1. 2014-15 (DHS-7) Kenya DHS: Kenya National Bureau of Statistics, Ministry of Health/Kenya,
 National AIDS Control Council/Kenya, Kenya Medical Research Institute, National Council for
 Population and Development/Kenya, and ICF International. 2015. Kenya Demographic and Health
 Survey 2014 [DATASETS]. Rockville, MD, USA: Kenya National Bureau of Statistics, Ministry
 of Health/Kenya, National AIDS Control Council/Kenya, Kenya Medical Research Institute,
 National Council for Population and Development/Kenya, and ICF International.

2. 2022 (DHS-8) Kenya DHS: KNBS and ICF. 2023. Kenya Demographic and Health Survey 2022:
 [DATASETS]. Nairobi, Kenya, and Rockville, Maryland, USA: KNBS and ICF. 
 <http://dhsprogram.com/pubs/pdf/FR339/FR339.pdf>. 2017

3.  The DHS Program Code Share Project, Code Library, DHS Program. DHS
    Program GitHub site. <https://github.com/DHSProgram>., in DHS
    Program GitHub site. 2022.

4. Rue, H., Martino, S. and Chopin, N., 2009. Approximate Bayesian inference
   for latent Gaussian models by using integrated nested Laplace approximations.
   Journal of the Royal Statistical Society Series B: Statistical Methodology, 71(2), pp.319-392.

5. Lindgren, F., Rue, H. and Lindström, J., 2011. An explicit link between Gaussian fields and Gaussian
   Markov random fields: the stochastic partial differential equation approach. Journal of the Royal
   Statistical Society Series B: Statistical Methodology, 73(4), pp.423-498.

6. WorldPop (www.worldpop.org - School of Geography and Environmental Science, University of Southampton;
 Department of Geography and Geosciences, University of Louisville; Departement de Geographie, Universite de Namur)
 and Center for International Earth Science Information Network (CIESIN), Columbia University (2018).
 Global High Resolution Population Denominators Project - Funded by The Bill and Melinda Gates Foundation (OPP1134076).
 https://dx.doi.org/10.5258/SOTON/WP00670

------------------------------------------------------------------------
