# Phylogenomics of Monitor Lizards

This repository holds associated files and scripts for our manuscript on the phylogenomics of *Varanus* lizards and the evolution of size disparity in a coevolutionary framework. Included are data files and code to replicate phylogenetic comparative analyses of body size evolution and functional diversity. In each folder is an RMarkdown document which will walk you through the desired analyses (*Trait* or *Spatial Evolution*). The easiest thing is to clone this repo and you should have everything you need categorized in the folders below. However, if you are looking for something that I've forgotten to put up, please email me at ian.brennan.anu.edu.au

Otherwise, this repository includes the following (organized under folders of the same names):

 ## Trait_Evo
 + **All_Size_Data.csv**: data file holding size (SVL--snout-vent-length; Tail), group (marsupial or varanoid), habitat, and current status (extinct or extant) information for all taxa included in this study.
 + **CoEvo_Walkthrough.RDS**: R data file which contains all the phylogenetic, size, and distributional data for the coevolutionary analyses of monitor lizards and marsupials in the *Monitor Phylogenomics* walkthrough.
 + **Goanna_Ancestors.gif**: animated gif of the spread and speciation of monitor lizards across the Australian continent.
 + **Goanna_Walkthrough.RDS**: R data file which contains all the phylogenetic, size, and distributional data for the evolutionary analyses of monitor lizards in the *Monitor Phylogenomics* walkthrough.
 + **GoannaMaps**: contemporary distribution maps of Australian monitors based on Atlas of Living Australia records. A subfolder holds the distribution maps of ancestral nodes.
 + **MonitorPhylogenomics** .html, .pdf, .Rmd: walkthrough files for the trait evolution of monitor lizards and marsupials. Exported as a webpage (html), pdf, and the base RMarkdown file.
 + **Varanus_Empirical_ModelObjects.RData** and **Varanus.Marsupial_Empirical_ModelObjects.RData**: if following along with the walkthrough, these files speed up the process by providing fitted model objects instead of requiring the user to run the entire model fitting exercises. 
