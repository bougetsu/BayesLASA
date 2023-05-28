### Pathology image data

This folder contains the input data for the analysis of tumor pathology images of lung cancer patients application case in the manuscript "Bayesian Landmark-based Shape Analysis of Tumor Pathology Images".

Processed polygonal chains of tumor region from 267 pathology images were provided in Rdata format names as `XXXX_outline.Rdata`. The raw tumor pathology images from the National Lung Screening Trial (NLST) are available at https://cdas.cancer.gov/datasets/nlst/

The Rdata contains

* `outline_polygon`: R list of processed polygonal chains of tumor region from pathology images. Each item in the list represents a tumor region in the sample and the two columns are x- and y- coordinates respectively.

* `outline_polygon`: R list of processed polygonal chains of tumor region from pathology images. Redundant vertices along the same straight segment were eliminated to only keep the start and end points. Each item in the list represents a tumor region in the sample and the two columns are x- and y- coordinates respectively.