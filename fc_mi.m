function signal_dataset = fc_mi(d, T, verbose)

%% Multi-Resolution-Analysis Deterministic Dendritic Cell Algorithm with Segmentation (MRA S-dDCA)
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
%% Function description.
% This function implements the feature selection approach used to perform
% signal categorization in the MRA S-dDCA model.

% Parameters:
% d: Dataset table, as provided by MRA_SdDCA function.
% T: Number of features to be selected for each signal category, as
%    provided by MRA_SdDCA function.
% verbose: Display algorithm relevant information.

if verbose
    disp("Obtain weight average of attack and normal datasets");
end
newid = 1:size(d.dataset_label,1);
d_temp = d;
for i=1:3
    d_temp(:, end) = [];   
end
dataset_normal = d_temp(d.dataset_label == 0,:);
dataset_attack = d_temp(d.dataset_label == 1,:);
mitable = table;
for i = 1 : size(d_temp,2)
    mitable.feature(i) = i;
    minormal = mi(table2array(dataset_normal(:,i)), d.cat);
    miattack = mi(table2array(dataset_attack(:,i)), d.cat); 
    mitable.normal(i) = minormal;
    mitable.attack(i) = miattack;
    diff = abs(minormal) - abs(miattack);
    mitable.diff(i) = diff;
    if diff > 0
        mitable.category(i) = 1;
        mitable.catname(i) = "SS";
    else
        mitable.category(i) = 0;  
        mitable.catname(i) = "DS";
    end  
    mitable.name(i) =  string(d.Properties.VariableNames{i});
end

%% Select best T separation features
mitable2 = table;
[~, Isafe] = maxk(table2array(mitable(:,4)),T);
[~, Idanger] = mink(table2array(mitable(:,4)),T);
contfsel = 1;
for i = 1 : T
    mitable2.feature(contfsel) = mitable.feature(Isafe(i));
    mitable2.normal(contfsel) = mitable.normal(Isafe(i));
    mitable2.attack(contfsel) = mitable.attack(Isafe(i));    
    mitable2.diff(contfsel) = mitable.diff(Isafe(i));
    mitable2.category(contfsel) = 1;
    mitable2.catname(contfsel) = "SS";    
    mitable2.name(contfsel) = mitable.name(Isafe(i));
    contfsel = contfsel + 1;  
end

for i = 1 : T
    mitable2.feature(contfsel) = mitable.feature(Idanger(i));
    mitable2.normal(contfsel) = mitable.normal(Idanger(i));
    mitable2.attack(contfsel) = mitable.attack(Idanger(i));    
    mitable2.diff(contfsel) = mitable.diff(Idanger(i));
    mitable2.category(contfsel) = 0;
    mitable2.catname(contfsel) = "DS";    
    mitable2.name(contfsel) = mitable.name(Idanger(i));
    contfsel = contfsel + 1;
end

mitable = mitable2;

%% Determine signal category (ss, pamp, ds)+++++++++++
if verbose
    disp("Determine signal category (ss, pamp, ds)+++++++++++");
end
featsel_1 = table;
featsel_1.feature = mitable.feature;
featsel_1.category = mitable.category;
featsel_1.name = mitable.name;
%% Preprocessing & initialization phase
if verbose
    disp("DCA: preprocessing & initialization phase");
end

% Signal category separation

ss_derived_f_1 = [];
ds_derived_f_1 = [];

for i = 1 : size(featsel_1, 1)
    if table2array(featsel_1(i,2)) == 0
        ds_derived_f_1 = [ds_derived_f_1; featsel_1(i,1)];
    elseif table2array(featsel_1(i,2)) == 1
        ss_derived_f_1 = [ss_derived_f_1; featsel_1(i,1)];
    end
end

%% Signal dataset generation
if verbose
    disp("Signal dataset generation");
end
% Danger signal generation

ss_derived_1 = zeros(size(d_temp, 1), 1);
for i = 1 : size(ss_derived_f_1,1)
    tmpel = table2array( d_temp(:,ss_derived_f_1.feature(i)) );
    ss_derived_1 = ss_derived_1 + tmpel;
end
ss_derived_1 = ss_derived_1 / size(ss_derived_f_1,1);

% SS and PAMP generation
ds_derived_1 = zeros(size(d_temp, 1), 1);
for i = 1 : size(ds_derived_f_1,1)
    tmpel = table2array( d_temp(:,ds_derived_f_1.feature(i)) );
    ds_derived_1 = ds_derived_1 + tmpel;
end
ds_derived_1 = ds_derived_1 / size(ds_derived_f_1,1);

% Final signal dataset to be used by DCA
signal_dataset = {newid', normalize(ds_derived_1, 'range'), normalize(ss_derived_1, 'range'), d.dataset_label, d.newid_addon};

end

