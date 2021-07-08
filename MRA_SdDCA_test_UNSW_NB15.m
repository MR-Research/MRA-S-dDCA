%% Multi-Resolution-Analysis Deterministic Dendritic Cell Algorithm with Segmentation (MRA S-dDCA)
% Testing script for the UNSW-NB15 dataset.
% Based on the work done by Greensmith in 2008 [1],
% adapted the segmentation concept proposed by Gu et al., 2009 [2].
% Developed by David Limon-Cantu, last modified May 2021.
% "Dendritic cells are immune sentinels"
% References:
% [1]J. Greensmith and U. Aickelin, “The Deterministic Dendritic Cell Algorithm,”
%    in Artificial Immune Systems, 2008, pp. 291–302.
% [2]F. Gu, J. Greensmith, and U. Aickelin, 
%    “Integrating Real-Time Analysis with the Dendritic Cell Algorithm through Segmentation,” 
%    in Proceedings of the 11th Annual Conference on Genetic and Evolutionary Computation, 
%    New York, NY, USA, 2009, pp. 1203–1210. doi: 10.1145/1569901.1570063.
% -------------------------------------------------------------------

% Load the dataset.
dataset_name = 'UNSW_NB15_testing-set.csv';
dataset = UNSWNB15Load(dataset_name);

% Prepare dataset.
newid_addon = strcat(string(table2array(dataset(:,3))), '_', string(table2array(dataset(:,4))), '_', string(table2array(dataset(:,5))));
dataset_label = dataset.label;
dataset_cat = dataset.attack_cat;
dataset = removevars(dataset, 'id');
dataset = removevars(dataset, 'proto');
dataset = removevars(dataset, 'service');
dataset = removevars(dataset, 'state');
dataset = removevars(dataset, 'attack_cat');
dataset = removevars(dataset, 'label');
dataset.cat = dataset_cat;
dataset.newid_addon = newid_addon;
dataset.dataset_label = dataset_label;

% MRA S-dDCA function call.
[confmat, accuracy, runtime] = MRA_SdDCA(dataset, "UNSW-NB15", 5, 2204, 128, 1, "sym2", false);