# Phylogenomics of Monitor Lizards

This repository holds associated files and scripts for our manuscript on the phylogenomics of *Varanus* lizards and the evolution of size disparity in a coevolutionary framework. Included are data files and code to replicate phylogenetic comparative analyses of body size evolution and functional diversity. In each folder is an RMarkdown document which will walk you through the desired analyses (*Trait* or *Spatial Evolution*). The easiest thing is to clone this repo and you should have everything you need categorized in the folders below. However, if you are looking for something that I've forgotten to put up, please email me at ian.brennan.anu.edu.au

Otherwise, this repository includes the following (organized under folders of the same names):

 ## Trait_Evo
 This folder is organized around the **Monitor Phylogenomics** walkthrough file (.Rmd) which includes code and explanation for plotting the phylogeny, manipulating data, and fitting models of trait evolution include coevolutionary models.
 + **All_Size_Data.csv**: data file holding size (SVL--snout-vent-length; Tail), group (marsupial or varanoid), habitat, and current status (extinct or extant) information for all taxa included in this study.
 + **CoEvo_Walkthrough.RDS**: R data file which contains all the phylogenetic, size, and distributional data for the coevolutionary analyses of monitor lizards and marsupials in the *Monitor Phylogenomics* walkthrough.
 + **Goanna_Ancestors.gif**: animated gif of the spread and speciation of monitor lizards across the Australian continent.
 + **Goanna_Walkthrough.RDS**: R data file which contains all the phylogenetic, size, and distributional data for the evolutionary analyses of monitor lizards in the *Monitor Phylogenomics* walkthrough.
 + **GoannaMaps**: contemporary distribution maps of Australian monitors based on Atlas of Living Australia records. A subfolder holds the distribution maps of ancestral nodes.
 + **MonitorPhylogenomics** .html, .pdf, .Rmd: walkthrough files for the trait evolution of monitor lizards and marsupials. Exported as a webpage (html), pdf, and the base RMarkdown file.
 + **Varanus_Empirical_ModelObjects.RData** and **Varanus.Marsupial_Empirical_ModelObjects.RData**: if following along with the walkthrough, these files speed up the process by providing fitted model objects instead of requiring the user to run the entire model fitting exercises. 

## Spatial_Evo
This folder is organized around the **Spatial Functional Diversity** walkthrough file (.Rmd) which includes code and explanation for visualizing and analyzing spatial and functional diversity of Australian monitors and faunivorous marsupials.
+ **Spatial_Walkthrough.RDS**, **Marsupial_Walkthrough.RDS** and **CoSpatial_Walkthrough.RDS**: R data files which contain all the phylogenetic and distributional data for the coevolutionary analyses of monitor lizards and marsupials in the *Spatial_Functional_Diversity* walkthrough.
+ **DispersalNullModel.R**: the code for the dispersal null model used to sample functional diversity across grid cells.
+  **SimulatedGoanna...**, **SimulatedMarsupial...**, and **SimulatedBoth...RDS**: R data files containing functional diversity data simulated using the *Dispersal Null Model* for each of the montitors, marsupials, and both together.
+ **Spatial_Functional_Diversity** .pdf, .Rmd: walkthrough files for the spatial evolution and functional diversity analyses. 

## Scripts_Files
This folder holds a number of functions and scripts used to analyze the data and that are sourced in the walkthrough files listed above.
+ **Calculate_AICs.R**: functions to quickly calculate AIC and AICc values to determine the best fitting models of trait evolution.
+ **CoEvo_Goanna_Distributions.csv** and **CoEvo_Marsupial_Distributions.csv**: hard copies of the distribution data for monitor lizards and marsupials used in this study.
+ **CreateGeoObject_fromSP.R**: designed to follow the *CreateGeoObject* function of RPANDA, this function creates geo_objects from explicit spatial data points and geometries. Better explained in the manuscript.
+ **make.OUwie.input_Script.R**: this helper function takes a csv file and makes an OUwie input object.
+ **plot.distmaps.R**: this function takes in spatial data as lat/lon and plots the distributions, while translating the data into spatial geometries.
+ **process.rase**: follows running *rase* and processes the output into a similar format as the empirical data, resulting in a similar output to *plot.distmaps*
+ **RPANDA_extras.R**: includes the additional evolutionary models described in the text. It also includes a number of RPANDA functions not included in the current CRAN version which are required to fit the models.
**search.surface.R**: is a relatively simple script to refit a model or models in parallel from different starting parameters to facilitate model fitting with complex models.

## Map_ShapeFiles
