# Multi-Resolution-Analysis Deterministic Dendritic Cell Algorithm with Segmentation (MRA S-dDCA).

Based on the work done by Greensmith in 2008 [1], adapted the segmentation concept proposed by Gu et al., 2009 [2]. Developed by David Limon-Cantu, last modified May 2021.

## Function Description.
This is a dDCA inspired intrusion detection model, developed as an Intrusion Detection System (IDS) approach.

## Parameters:
* **d**: High-level network feature dataset.
* **name**: Dataset name, used for storing results in file.
* **T**: Number of features to be selected for each signal category.
* **S**: Random number generator seed used to assign the migration threshold for the dendritic cell population.
* **m**: Segment size used by the S-dDCA.
* **p**: Dendritic cell population size.
* **w**: Wavelet used for the MODWT process.
* **verbose**: Display algorithm relevant information.

## Usage:
Provide the algorithm with your high-level network feature dataset "d" (such as the NSL-KDD, KDD99, UNSW-NB15).

> * You can get the UNSW-NB15 testing set here: https://research.unsw.edu.au/projects/unsw-nb15-dataset
> * You can get the NSL-KDD testing set here: https://www.unb.ca/cic/datasets/nsl.html

Provided features (dataset) "d" must be a table containing all (numeric) high-level network features available to the dataset. Non-real numeric values, such as categorical, logical, or string values are not allowed.
The last three features of the dataset  "d" must be as follows:
  > - "cat" must be the dataset attack labels (or classess)
  > - "newid_addon" must be a categorical variable that identifies the network flow instance (dataset row) using categorical or string values from the network dataset.
  > - "dataset_label" must be the binary dataset labels, where 0 denotes normal behavior, and 1 denotes (any) network anomaly.
  > 
> The resulting output is recorded in a CSV file following the name
convention: WAV_"name"-Y-M-D-h-m-s.ms.csv, where "name" is the name argument 
provided as function input, and the date and time of the performed test is used,
where Y is year, M is month, D is day, h is hour, m is minute, s is second, 
and ms is miliseconds.

## Additional functions:
* fc_mi.m implements the feature selection approach used to perform signal categorization in the MRA S-dDCA model. 
* mi.m and hist2.m determines the mutual information of two images or signals, and were developed in 2015 by Jose Delpiano, and modified by Andrew Hill in 2011. The corresponding license of these files is included in their respective files.

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see  [https://www.gnu.org/licenses/](https://www.gnu.org/licenses/).

## References:
[1]J. Greensmith and U. Aickelin, “The Deterministic Dendritic Cell Algorithm,”
   in Artificial Immune Systems, 2008, pp. 291–302.
[2]F. Gu, J. Greensmith, and U. Aickelin, 
   “Integrating Real-Time Analysis with the Dendritic Cell Algorithm through Segmentation,” 
   in Proceedings of the 11th Annual Conference on Genetic and Evolutionary Computation, 
   New York, NY, USA, 2009, pp. 1203–1210. doi: 10.1145/1569901.1570063.