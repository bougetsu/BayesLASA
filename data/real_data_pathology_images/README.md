## Pathology image dataset

This folder contains the tumor pathology image dataset taken from lung cancer patients for testing BayesLASA, as demonstrated in the manuscript "Bayesian Landmark-based Shape Analysis of Tumor Pathology Images".

Processed polygonal chains of tumor regions taken from 246 pathology images are provided in .Rdata format, as `XXXX_outline.Rdata`. The raw tumor pathology images from the National Lung Screening Trial (NLST) are available at https://cdas.cancer.gov/datasets/nlst/

Each .Rdata file contains:

* `outline_polygon`: An R List of processed polygonal chains of tumor regions from the pathology image. Each item in the list represents a tumor region in the sample and the two columns are x- and y-coordinates respectively.

* `outline_polygon`: An R List of processed polygonal chains of tumor regions from the pathology image. Redundant vertices along the same straight segment were eliminated to keep only the start and end points. Each item in the list represents a tumor region in the sample and the two columns are x- and y-coordinates respectively.
