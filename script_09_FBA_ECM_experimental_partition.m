clearvars; 
clc;

P3G_labeling_table = readtable('P3G_labeling_data.csv');
% Example MATLAB code to calculate the mean within a time window

% Define your time window (e.g., between 3 and 9 minutes)
time_window_start = 0; % in minutes
time_window_end = 21;   % in minutes

% Extract the relevant rows from the table
rows_in_window = (P3G_labeling_table.t_experimental >= time_window_start) & ...
                 (P3G_labeling_table.t_experimental <= time_window_end);

% Filter the 'mean' values based on the time window
filtered_mean_values = P3G_labeling_table.mean(rows_in_window);

% Calculate the mean of the filtered values
median_within_window = median(filtered_mean_values);

% Display the result
fprintf('The median percent of 13C in the plateau between %.1f and %.1f minutes is %.4f.\n', ...
        time_window_start, time_window_end, median_within_window);

percent_13C_P3G_plateau = median_within_window;
C13max = (1/6)*100;
percent_P3G_from_rubisco =  (percent_13C_P3G_plateau/C13max)*100;
percent_P3G_from_pgk = 100 - percent_P3G_from_rubisco;

%% defining the experimental data to be considered
% labeling = 1; cycle = 2
experiment = 1;
%% defining the saturation state of the enzymes
% KM=S sst=1; KM=S/10 sst=2
saturation_state = 1; 
%% assigning state
if experiment == 1
   qAc = -0.52;
   experiment_label={'labeling'};
else
   filename = 'reconciled_rates_cycle.csv';
   dataTable = readtable(filename);
   metabolite_names = dataTable{:, 1};
   associated_numbers = dataTable{:, 2};
   target_metabolite = 'acetate';
   metabolite_index = find(strcmp(target_metabolite, metabolite_names));
   qAc = associated_numbers(metabolite_index);
   experiment_label={'cycle'};
end
    
if saturation_state == 1
   timesKM = 1;
   saturation_label = {'_KM_'};
else
   timesKM = 10;
   saturation_label = {'_sat_'};
end

% Define symbolic variables
syms mol_P3G_NTS mol_P3G_EMP percent_EMP

% Define the system of equations
eq1 = mol_P3G_NTS / mol_P3G_EMP == percent_P3G_from_rubisco / (100-percent_P3G_from_rubisco); % Equation I
eq2 = (percent_EMP / (100-percent_EMP)) == (mol_P3G_EMP * (0.5 * 6)) / (mol_P3G_NTS * 2.5); % Equation II
eq3 = mol_P3G_NTS + 2*mol_P3G_EMP == 0.5*((-1*qAc) + mol_P3G_NTS + mol_P3G_EMP); % Equation III

% Solve the system of equations
[sol_mol_P3G_NTS, sol_mol_P3G_EMP, sol_percent_EMP] = ...
solve([eq1, eq2, eq3], [mol_P3G_NTS, mol_P3G_EMP, percent_EMP]);

% Display the results
fprintf('mol of P3G formed in NTS: %.4f\n', double(sol_mol_P3G_NTS));
fprintf('mol of P3G formed in EMP: %.4f\n', double(sol_mol_P3G_EMP));
fprintf('percentage of glucose oxidized by NTS: %.4f\n', 100-double(sol_percent_EMP));
fprintf('percentage of glucose oxidized by EMP: %.4f\n', double(sol_percent_EMP));
y_horizontal_line = 100-double(sol_percent_EMP);
% line_plateau = y_horizontal_line;
mol_P3G_PGK = double(sol_mol_P3G_EMP);
mol_P3G_Rubisco = double(sol_mol_P3G_NTS);

disp ('-----------------------------')
disp ('The glycolytic partition is:')
percent_NTS = 100-double(sol_percent_EMP)
percent_EMP = double(sol_percent_EMP)

Table_output=zeros(1,3);
Table_output(1,1) = percent_EMP;

%% Flux Balance Analysis

order = {'GlycP';'GPM';'PGI';'PPiPFK';'FBPald';'TPI';'GAPDH';'PGK';'TktA';'SBPald';'SBPase';'TktB';'RibE';'RibI';'RbuK';'RbuCO';'PGM';'ENO';'PYK';'PDH';'ACS';'Thio';'AAR'}; % the desired order

load Accumulibacter_anaerobic;
model = changeRxnBounds(model,{'AMPPPT','ADPPPT','PYK','PDH','PHBsyn','SBPase','RbuK','RbuCO','ACS'},0,'l');
model = changeObjective(model,'EX_PHB');%selecting Objective
model = changeRxnBounds(model,'EX_Ace',qAc,'l');
model = changeRxnBounds(model,'EX_Mal4',-1000,'l');
model = changeRxnBounds(model,'EX_Mal3',1000,'u');
model = changeRxnBounds(model,'EX_PHB',1000,'u');
model = changeRxnBounds(model,'RbuCO',double((sol_mol_P3G_NTS))/2,'b');
FBAsolution = optimizeCbModel(model,'max');

qGlyc = FBAsolution.x(findRxnIDs(model,'EX_Mal4'));
Table_output(1,2) = qGlyc;
Table_output(1,3) = qAc;

for i=1:length(order)
    normalized_fluxes(i) = FBAsolution.x( findRxnIDs(model,order(i)) )/(-1*qAc);      
        if normalized_fluxes(i)<0.000001
            normalized_fluxes(i) = 0;
        end
end

ActiveRxnsFormula = printRxnFormula(model,order,false); %clc;

% Define the string to be replaced and the replacement string
stringToReplace = {'->'};
replacementString = {'<=>'};

% Iterate through the array of strings and replace the target string
for i = 1:numel(ActiveRxnsFormula)
    ActiveRxnsFormula{i} = strrep(ActiveRxnsFormula{i}, stringToReplace, replacementString);
end

% Define the string to be replaced and the replacement string
stringToReplace = {'  '};
replacementString = {' '};

% Iterate through the array of strings and replace the target string
for i = 1:numel(ActiveRxnsFormula)
    ActiveRxnsFormula{i} = strrep(ActiveRxnsFormula{i}, stringToReplace, replacementString);
end

%% eliminating reactions that are not analyzed in eQuilibrator

list_zero_flux = {};

for i=1:length(model.rxns)
    if FBAsolution.x(i)==0
       list_zero_flux = [list_zero_flux, model.rxns(i)];
    end
end

exchange_rxns = findExcRxns(model,1,1);

list_exchange = {};

for i=1:length(model.rxns)
    if exchange_rxns(i)==1
       list_exchange = [list_exchange, model.rxns(i)];
    end
end

rxns_to_remove = [list_exchange,list_zero_flux,'PPiase','AMPPPT','ADPPPT','PHBsyn'];

model = removeRxns(model,rxns_to_remove);

%% getting data obtained with eQuilibrator

%%% Get the metabolite concentrations calculated by MDF analysis
file_name = 'compound_concentrations.csv';
delimiter = "	";
concentrations_table = readtable(file_name, 'Delimiter', delimiter);
%re-ordering the metabolites in the same order of the stoichiometric model
[~, idx] = ismember(concentrations_table.Compound, model.mets);
[~, sortorder] = sort(idx);
concentrations_table = concentrations_table(sortorder,:);
list_mets = table2array(concentrations_table(:,1));
concs_MDF = table2array(concentrations_table(:,2));

%%% Get the standard deltaG of the reactions participating in the network
file_name = 'standard_dGs.csv';
delimiter = "	";
dGs_table = readtable(file_name, 'Delimiter', delimiter);
%re-ordering the reactions in the same order of the stoichiometric model
[~, idx] = ismember(dGs_table.reactions, model.rxns);
[~, sortorder] = sort(idx);
dGs_table = dGs_table(sortorder,:);
list_rxns = table2array(dGs_table(:,1));
list_dGs = table2array(dGs_table(:,2));
R = 0.008314; T = 20 + 273.15;
%calculating equilibrium constants from standard dGs 
list_Keqs = exp((list_dGs.*-1)/(R*T));

%%% Creating the matrices to be filled
KMS=zeros(length(model.mets),length(model.rxns));
stoich_coef_substrates=zeros(length(model.mets),length(model.rxns));
KMP=zeros(length(model.mets),length(model.rxns));
stoich_coef_products=zeros(length(model.mets),length(model.rxns));
kcatR = zeros(length(model.rxns),1);
kcatF = zeros(length(model.rxns),1);
weights = zeros(length(model.rxns),1);

% Read the table from the file
protein_table = readtable('enzymes_kcats_MWs_medians.csv');

% assigning the kinetic parameters

    for r=1:length(model.rxns)

        for m=1:length(model.mets)

            if model.S(m,r)<0

                    KMS(m,r)=(concs_MDF(m)*1000/timesKM); %converting to mM
                    stoich_coef_substrates(m,r)=abs(model.S(m,r));

                if strcmp(model.metNames{m}, 'H2O')
                    KMS(m, r) = 1;
                    stoich_coef_substrates(m, r) = abs(model.S(m, r));
                end

            end

            if model.S(m,r)>0
                KMP(m,r)=(concs_MDF(m)*1000); %converting to mM
                stoich_coef_products(m,r)=abs(model.S(m,r));

                if strcmp(model.metNames{m}, 'H2O')
                    KMP(m, r) = 1;
                    stoich_coef_products(m, r) = abs(model.S(m, r));
                end
            end

        end

      %for each reaction, get the KM values from the matrix KMS 
      A = KMS(:,r);
      valuesKMS = A(A ~= 0);
      C = stoich_coef_substrates(:,r);
      valuescoef_substrates = C(C ~= 0);
      %Each KMS value is powered to its corresponding stoichiometric coefficient
      adjusted_substrates = valuesKMS .^ valuescoef_substrates;

      B = KMP(:,r);
      valuesKMP = B(B ~= 0);
      order2 = stoich_coef_products(:,r);
      valuescoef_products = order2(order2 ~= 0);
      %Each KMP value is powered to its corresponding stoichiometric coefficient
      adjusted_products = valuesKMP .^ valuescoef_products;

        %Specify the reaction name you're interested in
        target_reaction_name = model.rxns(r);  % Replace with the desired reaction name

        %Find the row index for the specified reaction name
        row_index = find(strcmp(protein_table.Reaction_Name, target_reaction_name));
        weights(r) = protein_table.MW(row_index);

        if ismember(target_reaction_name, {'TktA','TktB'})
            kcatR(r) = protein_table.kcat(row_index);
            kcatF(r) = kcatR(r) * (list_Keqs(r) * prod(adjusted_substrates)) / prod(adjusted_products);
        else
            % Get the corresponding kcat values
            kcatF(r) = protein_table.kcat(row_index);
            % kcat reverse are calculated using the Haldane relationship
            kcatR(r) = (kcatF(r) * prod(adjusted_products)) / (list_Keqs(r) * prod(adjusted_substrates));
        end

    end
    
%exporting the kinetic parameters in the required format
combined_matrix = KMS+KMP;
[metabolites,reactions]=size(combined_matrix);
chain = [];
chain2 = [];
KMs = cell(reactions,1);

for g=1:reactions
    for u=1:metabolites
        if combined_matrix(u,g)~0;
           chain = strcat(model.mets(u),':',num2str(combined_matrix(u,g)));
           chain2 = strcat(chain2,{' '},chain);
        end          
    end
        KMs(g) = chain2;
        chain2 = [];
end

all_KMS = cell2table(KMs);
T1 = table(model.rxns,kcatF,kcatR);
T4 = table(weights);
T1 = [T1, all_KMS, T4];
T1.Properties.VariableNames = {'Reaction Name' 'kcatf (1/s)' 'kcatr (1/s)' 'KM (mM)' 'MWe(Da)'};
[~, idx] = ismember(T1.("Reaction Name"), order);
[~, sortorder] = sort(idx);
T1 = T1(sortorder,:);    
file_name = strcat('consistent_kinetic_parameters',string(saturation_label),'.csv');
writetable(T1,file_name);
[~, idx] = ismember(protein_table.("Reaction_Name"), order);
[~, sortorder] = sort(idx);
protein_table = protein_table(sortorder,:);
T2 = table(ActiveRxnsFormula);
T2.Properties.VariableNames = {'Reaction Formula'}; 
T3 = table(normalized_fluxes');
T3.Properties.VariableNames = {'Relative Flux'}; 
T5 = [T2,T3,T1];
T5.Properties.VariableNames{'kcatf (1/s)'} = 'kcrf(1/s)';
T5.Properties.VariableNames{'kcatr (1/s)'} = 'kcrr(1/s)';
T5.Properties.VariableNames{'KM (mM)'} = 'kM(mM)';
pathway_name = strcat('NTS_and_EMP_meta',string(saturation_label),string(experiment_label),'.csv');      
% Remove rows where "Relative Flux" is zero
T5 = T5(T5.('Relative Flux') ~= 0, :);
% Save the updated table as a CSV file
writetable(T5, string(pathway_name));
% T6 = array2table(Table_output,'VariableNames',{'Percent NTS' 'qGlycogen' 'qAcetate'});
% writetable(T6,'qrates_experimental_partition.csv');