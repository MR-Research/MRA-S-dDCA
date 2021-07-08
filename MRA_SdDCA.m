function [confmat, a4_accuracy, timing] = MRA_SdDCA(d, name, T, S, m, p, w, verbose)

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
% This is a dDCA inspired intrusion detection model, developed as an
% Intrusion Detection System (IDS) approach.

% Parameters:
% d: High-level network feature dataset.
% name: Dataset name, used for storing results in file.
% T: Number of features to be selected for each signal category.
% S: Random number generator seed used to assign the migration.
% threshold for the dendritic cell population.
% m: Segment size used by the S-dDCA.
% p: Dendritic cell population size.
% w: Wavelet used for the MODWT process.
% verbose: Display algorithm relevant information.

% Usage:
% Provide the algorithm with your high-level network feature dataset "d" 
% (such as the NSL-KDD, KDD99, UNSW-NB15).
% Provided features (dataset) "d" must be a table containing all (numeric) 
% high-level network features available to the dataset. Non-real numeric values,
% such as categorical, logical, or string values are not allowed*.
% * The last three features of the dataset  "d" must be as follows:
%   - "cat" must be the dataset attack labels (or classess)
%   - "newid_addon" must be a categorical variable that identifies the
%     network flow instance (dataset row) using categorical or string values
%     from the network dataset.
%   - "dataset_label" must be the binary dataset labels, where 0 denotes
%     normal behavior, and 1 denotes (any) network anomaly.
% The resulting output is recorded in a CSV file following the name
% convention: WAV_"name"-Y-M-D-h-m-s.ms.csv, where "name" is the name argument 
% provided as function input, and the date and time of the performed test is used,
% where Y is year, M is month, D is day, h is hour, m is minute, s is second, 
% and ms is miliseconds.

% Additional functions:
% * fc_mi.m implements the feature selection approach used to perform signal 
%   categorization in the MRA S-dDCA model. 
% * mi.m and hist2.m determines the mutual information of two images or signals,
%   and were developed in 2015 by Jose Delpiano, and modified
%   by Andrew Hill in 2011. The corresponding license of these files is
%   included in their respective files.

arguments
    d (:,:) table
    name (1,1) string
    T (1,1) double
    S (1,1) double
    m (1,1) double
    p (1,1) double
    w (1,1) string
    verbose (1,1) logical = true
end

tmpd = clock;
rng(S);
filen = strcat("WAV_", name, '-', string(tmpd(1)),'-',string(tmpd(2)),'-',string(tmpd(3)),'-',string(tmpd(4)),'-',string(tmpd(5)),'-',string(tmpd(6)));
writecell({'Dataset' 'Wavelet' 'Wavelet level' 'DC Pop Size' 'Segment Size' 'Selected Features each Category' 'DCA version' 'Observations' 'Runtime' 'TP' 'TN' 'FP' 'FN' 'Precision' 'Sensitivity' 'Specificity' 'Accuracy'}, strcat(filen,'.csv'));
tic
% Feature-class mutual information maximization-based method for feature
% categorization and selection.
signal_dataset = fc_mi(d, T, verbose);

err = 0;
try
    wavlvl = p;
    m_og = m;
    antigencat = unique(unique(signal_dataset{1}));
    tp = 0;
    tn = 0;
    fp = 0;
    fn = 0;
    khtemp = zeros(1,p);
    csmtemp = zeros(1,p);
    alphacount = zeros(p,m);
    alpha = zeros(1,size(antigencat,1));
    kalpha = zeros(1,size(antigencat,1));      
    mt = rand(1,p);  
    if verbose
        disp("dDCA start...");  
    end
    list = zeros(1,m);
    signal_block = [signal_dataset{2}(1:m) signal_dataset{3}(1:m)];
    % MODWT of the first signal segment, pefrormed for both signal categories, namely DS and SS.
    wav_signal_ds = modwt(signal_block(:,1), w, wavlvl);
    wav_signal_ds = wav_signal_ds(1:end-1,:);
    wav_signal_ds = normalize(wav_signal_ds, 'range');
    wav_signal_ss = modwt(signal_block(:,2), w, wavlvl);
    wav_signal_ss = wav_signal_ss(1:end-1,:);
    wav_signal_ss = normalize(wav_signal_ss, 'range');
    % Signal energy calculation for current signal segment and decomposition
    % levels.
    ds_energy = (sum(wav_signal_ds.^2,2))';
    ss_energy = (sum(wav_signal_ss.^2,2))';
    energies_ds = zeros(size(signal_dataset{1},1),p);
    energies_ss = zeros(size(signal_dataset{1},1),p);
    classification = [signal_dataset{4}(:) zeros(size(signal_dataset{4},1),1)];
    last = -1;
    for i = 1 : size(signal_dataset{1},1)
        if (mod(i,1000) == 0 || i == size(signal_dataset{1},1)) && verbose                    
            disp(strcat("MRA S-dDCA iteration:", int2str(i), " of ", int2str(size(signal_dataset{1},1))));
        end
        % Segment processing.
        if (mod(i,m) ~= 0 && i <= size(signal_dataset{1},1)) || (i/m == floor(i/m))            
            if last == -1
                if mod(i,m) == 0
                    ind = m;
                else
                    ind = mod(i,m);
                end   
            else
                last = last+1;
                ind = last;
            end
            % dDCA detection phase.
            si = [wav_signal_ds(1:p,ind) wav_signal_ss(1:p,ind)];
            energies_ds(i,1:size(ds_energy,2)) = ds_energy;
            energies_ss(i,1:size(ss_energy,2)) = ss_energy;
            migrated = csmtemp <= mt;
            khtemp(migrated) = khtemp(migrated)+(si(migrated,1) - (2 *si(migrated,2)))';            
            alphacount(migrated,ind) = alphacount(migrated,ind) + 1;                    
            csmtemp(migrated) = csmtemp(migrated) + (si(migrated,1) + si(migrated,2))';
            kelem = find(csmtemp > mt);
            if length(kelem) >= 1
                alpha(i) = alpha(i) + sum(alphacount(kelem,ind));   
                kalpha(i) = kalpha(i) + sum(khtemp(kelem));
            end           
            khtemp(kelem) = 0;
            alphacount(kelem,:) = 0;
            csmtemp(kelem) = 0;
            if mod(i,m) ~= 0
                list(mod(i,m)) = i;
            else
                list(m) = i;
            end
        end
        % Segment processing finished.
        if (mod(i,m) == 0 && i < size(signal_dataset{1},1) && last == -1)  
            if i+m <= size(signal_dataset{1},1)
                signal_block = [signal_dataset{2}(i+1:i+m) signal_dataset{3}(i+1:i+m)];
            else
                signal_block = [signal_dataset{2}(i+1:size(signal_dataset{1},1)) signal_dataset{3}(i+1:size(signal_dataset{1},1))];
                wavlvl = floor(log2(size(signal_dataset{1},1)-(i +1)));
                m = size(signal_dataset{1},1)-(i);
                last = 0;
                if p == 1
                    mt = mt(1);
                else
                    if p > wavlvl
                        mt = mt(1:wavlvl);
                    else
                        mt = mt(1:p);
                    end
                end
                if wavlvl < p                                
                    p = wavlvl;
                end
                khtemp = khtemp(1:p);
                csmtemp = csmtemp(1:p);
                alphacount = alphacount(1:p,:);
            end
            % Next segment MODWT signal calculation.
            wav_signal_ds = modwt(signal_block(:,1), w, wavlvl);
            wav_signal_ds = wav_signal_ds(1:end-1,:);
            wav_signal_ds = normalize(wav_signal_ds, 'range');
            wav_signal_ss = modwt(signal_block(:,2), w, wavlvl);
            wav_signal_ss = wav_signal_ss(1:end-1,:);
            wav_signal_ss = normalize(wav_signal_ss, 'range');
            ds_energy = (sum(wav_signal_ds.^2,2))';
            ss_energy = (sum(wav_signal_ss.^2,2))';        
            alphacount = zeros(p,m);
            list = zeros(1,m);
        end
    end
    % DT classification
    if verbose
        disp("DCA DT classification model training...");
    end
    dcares = table;
    tmpka = (kalpha./alpha)';
    nanfilt = isnan(tmpka);
    tmpcat = d.newid_addon;
    dcares.ka = tmpka(~nanfilt);
    dcares.cat = tmpcat(~nanfilt);
    dcares.dsenergy = mean(energies_ds(~nanfilt,:),2);
    dcares.ssenergy = mean(energies_ss(~nanfilt,:),2);
    label_filt = d.dataset_label(~nanfilt);
    DTModel = fitctree(dcares, label_filt);
    if verbose
        disp("DCA DT classification model training finished...");
        disp("DCA testing, running trained DT model...");
    end
    
    prd = table;
    prd.ka = (kalpha./alpha)';
    prd.cat = d.newid_addon;
    prd.dsenergy = mean(energies_ds,2);
    prd.ssenergy = mean(energies_ss,2);
    label = predict(DTModel, prd);
    classification(:,2) = label; 
    
    if verbose
        disp("DT model testing finished...");
    end

    timing = toc;
    timing

    %% Metric generation

    for i = 1 : size(signal_dataset{1},1)
        if classification(i,1) == 0 && classification(i,1) == classification(i,2)
            tn = tn + 1;
        end

        if classification(i,1) == 1 && classification(i,1) == classification(i,2)
            tp = tp + 1;
        end    

        if classification(i,1) == 0 && classification(i,1) ~= classification(i,2)
            fn = fn + 1;
        end    

        if classification(i,1) == 1 && classification(i,1) ~= classification(i,2)
            fp = fp + 1;
        end    
    end
    
    if verbose
        disp("dDCA end...");
    end
    confmat = [tp tn fp fn];
    a4_accuracy = (tp+tn)/(tp+tn+fp+fn);
    a1_precision = (tp)/(tp+fp);
    a2_sensitivity = (tp)/(tp+fn);
    a3_specificity = (tn)/(tn+fp);
    
catch exception               
    err = 1;
    timing = toc;   
    if verbose
        timing
    end
    rethrow(exception)
end
% Write results to file            
if err == 0
    writecell({name w p p m_og T "1.0" size(signal_dataset{1},1) timing tp tn fp fn a1_precision a2_sensitivity a3_specificity a4_accuracy}, strcat(filen,'.csv'), 'WriteMode','append');
else
    writecell({name w p p m_og T "1.0" size(signal_dataset{1},1) timing "N/A" "N/A" "N/A" "N/A" "Error" "Error" "Error" "Error"}, strcat(filen,'.csv'), 'WriteMode','append');                
end 

end

