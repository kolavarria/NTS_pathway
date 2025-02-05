clearvars;
clc;

P3G_labeling_table = readtable('P3G_labeling_data.csv');

% Define your time window
time_window_start = 0; % in minutes
time_window_end = 21;   % in minutes

% Extract the relevant rows from the table
rows_in_window = (P3G_labeling_table.t_experimental >= time_window_start) & ...
                 (P3G_labeling_table.t_experimental <= time_window_end);

% Filter the 'mean' values based on the time window
filtered_values = P3G_labeling_table.mean(rows_in_window);

% Calculate the mean of the filtered values
median_within_window = median(filtered_values);

% Display the result
fprintf('The median percent of 13C in the plateau between %.1f and %.1f minutes is %.4f.\n', ...
        time_window_start, time_window_end, median_within_window);

percent_13C_P3G_plateau = median_within_window;
qAc=-0.52;
C13max = (1/6)*100;
percent_P3G_from_rubisco =  (percent_13C_P3G_plateau/C13max)*100;
percent_P3G_from_pgk = 100 - percent_P3G_from_rubisco;

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

%% Flux Balance Analysis
load Accumulibacter_anaerobic;
model = changeRxnBounds(model,{'ADPPPT','AMPPPT','PYK','PDH','PHBsyn','SBPase','RbuK','RbuCO','ACS'},0,'l');
model = changeObjective(model,'EX_PHB');%selecting Objective
model = changeRxnBounds(model,'EX_Ace',qAc,'l');
model = changeRxnBounds(model,'EX_Mal4',-1000,'l');
model = changeRxnBounds(model,'EX_Mal3',1000,'u');
model = changeRxnBounds(model,'EX_PHB',1000,'u');
model = changeRxnBounds(model,'RbuCO',double((sol_mol_P3G_NTS))/2,'b');
FBAsolution = optimizeCbModel(model,'max');
% disp ('----- Balanced q-rates associated to the glycolytic routes --------------')
% printFluxVector(model, FBAsolution.x,1,1);
ActiveRxnsFormula = printRxnFormula(model,model.rxns,false); %clc;
T1 = table(model.rxns,ActiveRxnsFormula,FBAsolution.x,'VariableNames',{'Reaction Name' 'Reaction Formula' 'Flux'});
% disp(T1)
%writetable(T1,'glycolytic_flux_distribution.txt','Delimiter','tab');
writetable(T1,'glycolytic_fluxes_labeling.csv');

%glycolytic q-rates
qGluc = FBAsolution.x(findRxnIDs(model,'EX_Mal4'));
qPi  = FBAsolution.x(findRxnIDs(model,'EX_Pi'));
qAce  = FBAsolution.x(findRxnIDs(model,'EX_Ace'));
qCO2  = FBAsolution.x(findRxnIDs(model,'EX_CO2'));
qHB  = FBAsolution.x(findRxnIDs(model,'EX_PHB'));
T2    = table(qGluc',qAce',qPi',qCO2',qHB','VariableNames',{'qGluc' 'qAce' 'qPi' 'qCO2' 'qHB'});
%writetable(T2,'q_rates_during_labeling.txt','Delimiter','tab');
writetable(T2,'q_rates_labeling.csv');

%% Generating the file containing the reactions carrying fluxes

order = {'ADPPPT';'AMPPPT';'GlycP';'GPM';'PGI';'PPiPFK';'FBPald';'TPI';'GAPDH';'PGK';'TktA';'SBPald';'SBPase';'TktB';'RibE';'RibI';'RbuK';'RbuCO';'PGM';'ENO';'PYK';'PDH';'ACS';'PPiase';'Thio';'AAR';'PHBsyn'}; % the desired order

for i=1:length(order)
    PHB_normalized_fluxes(i) = FBAsolution.x(findRxnIDs(model,order(i)))/FBAsolution.x(findRxnIDs(model,'EX_PHB'));
    glucose_normalized_fluxes(i) = FBAsolution.x(findRxnIDs(model,order(i)))/FBAsolution.x(findRxnIDs(model,'EX_Mal3'));
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

T3 = table(ActiveRxnsFormula,PHB_normalized_fluxes',order,'VariableNames',{'Reaction Formula' 'Relative Flux' 'Reaction Name'});
%writetable(T3,'PHB_normalized_glycolytic_fluxes.txt','Delimiter','tab');
writetable(T3,'PHB_normalized_fluxes_labeling.csv');

T4 = table(ActiveRxnsFormula,glucose_normalized_fluxes',order,'VariableNames',{'Reaction Formula' 'Relative Flux' 'Reaction Name'});
%writetable(T4,'glucose_normalized_glycolytic_fluxes.txt','Delimiter','tab');
writetable(T4,'glucose_normalized_fluxes_labeling.csv');

writetable(T4, 'NTS_and_EMP_FBA_solution_labeling.csv');

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

%%
rxns_to_remove = [list_exchange,list_zero_flux,'ADPPPT','AMPPPT','PHBsyn'];
%rxns_to_remove = [list_exchange,list_zero_flux];

% Identify rows to remove
rows_to_remove = ismember(T4.("Reaction Name"), rxns_to_remove);

% Remove rows from the table
T4(rows_to_remove, :) = [];

[~, idx] = ismember(T4.("Reaction Name"), order);
[~, sortorder] = sort(idx);
T4 = T4(sortorder,:);

writetable(T4, 'NTS_and_EMP_for_MDF.csv');

%% Visual flux partition

% Read the input data
P3G_labeling_table = readtable('P3G_labeling_data.csv');
mean_experimental = P3G_labeling_table.mean; % Renamed from "R_experimental"
t_experimental = P3G_labeling_table.t_experimental;
SD_experimental = P3G_labeling_table.SD;    % Renamed from "SD_R_experimental"

% Preallocate storage vectors for solutions
percent_NTS_vector = nan(length(mean_experimental), 1);

for i = 1:length(mean_experimental)
    
    percent_13C_P3G_plateau = mean_experimental(i);
    percent_P3G_from_rubisco = (percent_13C_P3G_plateau / C13max) * 100;
    percent_P3G_from_pgk = 100 - percent_P3G_from_rubisco;

    % Define symbolic variables
    syms mol_P3G_NTS mol_P3G_EMP percent_EMP

    % Define the system of equations
    eq1 = mol_P3G_NTS / mol_P3G_EMP == percent_P3G_from_rubisco / (100-percent_P3G_from_rubisco); % Equation I
    eq2 = (percent_EMP / (100-percent_EMP)) == (mol_P3G_EMP * (0.5 * 6)) / (mol_P3G_NTS * 2.5); % Equation II
    eq3 = mol_P3G_NTS + 2*mol_P3G_EMP == 0.5*((-1*qAc) + mol_P3G_NTS + mol_P3G_EMP); % Equation III

    % Solve the system of equations
    [sol_mol_P3G_NTS, sol_mol_P3G_EMP, sol_percent_EMP] = ...
        solve([eq1, eq2, eq3], [mol_P3G_NTS, mol_P3G_EMP, percent_EMP]);    
    
    % Check if solutions exist and convert them to numeric
    if ~isempty(sol_mol_P3G_NTS)
        % Convert solutions to numeric
        mol_P3G_NTS = double(sol_mol_P3G_NTS);
        percent_NTS = 100-double(sol_percent_EMP);
        
        % Handle multiple solutions: Choose the first valid solution
        if numel(mol_P3G_NTS) > 1
            mol_P3G_NTS = mol_P3G_NTS(1);
        end
        if numel(percent_NTS) > 1
            percent_NTS = percent_NTS(1);
        end
        
        % Save valid solutions in vectors
        mol_P3G_NTS_vector(i) = mol_P3G_NTS;
        percent_NTS_vector(i) = percent_NTS;
        
    else
        % If no solutions, assign NaN to indicate missing value
        mol_P3G_NTS_vector(i) = 0;
        percent_NTS_vector(i) = 0;
    end
end

%%%%%%%%%%%%%%%plus SD

% Preallocate storage vectors for solutions
percent_NTS_vector_plus = nan(length(mean_experimental), 1);

for i = 1:length(mean_experimental)
    percent_13C_P3G_plateau = mean_experimental(i) + SD_experimental(i);
    percent_P3G_from_rubisco = (percent_13C_P3G_plateau / C13max) * 100;
    percent_P3G_from_pgk = 100 - percent_P3G_from_rubisco;

 % Define symbolic variables
    syms mol_P3G_NTS mol_P3G_EMP percent_EMP

    % Define the system of equations
    eq1 = mol_P3G_NTS / mol_P3G_EMP == percent_P3G_from_rubisco / (100-percent_P3G_from_rubisco); % Equation I
    eq2 = (percent_EMP / (100-percent_EMP)) == (mol_P3G_EMP * (0.5 * 6)) / (mol_P3G_NTS * 2.5); % Equation II
    eq3 = mol_P3G_NTS + 2*mol_P3G_EMP == 0.5*((-1*qAc) + mol_P3G_NTS + mol_P3G_EMP); % Equation III

    % Solve the system of equations
    [sol_mol_P3G_NTS, sol_mol_P3G_EMP, sol_percent_EMP] = ...
        solve([eq1, eq2, eq3], [mol_P3G_NTS, mol_P3G_EMP, percent_EMP]);
    
    % Check if solutions exist and convert them to numeric
    if ~isempty(sol_mol_P3G_NTS)
        % Convert solutions to numeric
        mol_P3G_NTS = double(sol_mol_P3G_NTS);
        percent_NTS = 100-double(sol_percent_EMP);
        
        % Handle multiple solutions: Choose the first valid solution
        if numel(mol_P3G_NTS) > 1
            mol_P3G_NTS = mol_P3G_NTS(1);
        end
        if numel(percent_NTS) > 1
            percent_NTS = percent_NTS(1);
        end
        
        % Save valid solutions in vectors
        mol_P3G_NTS_vector(i) = mol_P3G_NTS;
        percent_NTS_vector_plus(i) = percent_NTS;
        
    else
        % If no solutions, assign NaN to indicate missing value
        mol_P3G_NTS_vector(i) = 0;
        percent_NTS_vector_plus(i) = 0;
    end
end

%%%%%%%%%%%%%%minus SD

% Preallocate storage vectors for solutions
percent_NTS_vector_minus = nan(length(mean_experimental), 1);

for i = 1:length(mean_experimental)
    percent_13C_P3G_plateau = mean_experimental(i) - SD_experimental(i);
    percent_P3G_from_rubisco = (percent_13C_P3G_plateau / C13max) * 100;
    percent_P3G_from_pgk = 100 - percent_P3G_from_rubisco;


 % Define symbolic variables
    syms mol_P3G_NTS mol_P3G_EMP percent_EMP

    % Define the system of equations
    eq1 = mol_P3G_NTS / mol_P3G_EMP == percent_P3G_from_rubisco / (100-percent_P3G_from_rubisco); % Equation I
    eq2 = (percent_EMP / (100-percent_EMP)) == (mol_P3G_EMP * (0.5 * 6)) / (mol_P3G_NTS * 2.5); % Equation II
    eq3 = mol_P3G_NTS + 2*mol_P3G_EMP == 0.5*((-1*qAc) + mol_P3G_NTS + mol_P3G_EMP); % Equation III

    % Solve the system of equations
    [sol_mol_P3G_NTS, sol_mol_P3G_EMP, sol_percent_EMP] = ...
        solve([eq1, eq2, eq3], [mol_P3G_NTS, mol_P3G_EMP, percent_EMP]);

    % Check if solutions exist and convert them to numeric
    if ~isempty(sol_mol_P3G_NTS)
        % Convert solutions to numeric
        mol_P3G_NTS = double(sol_mol_P3G_NTS);
        percent_NTS = 100-double(sol_percent_EMP);
        
        % Handle multiple solutions: Choose the first valid solution
        if numel(mol_P3G_NTS) > 1
            mol_P3G_NTS = mol_P3G_NTS(1);
        end
        if numel(percent_NTS) > 1
            percent_NTS = percent_NTS(1);
        end
        
        % Save valid solutions in vectors
        mol_P3G_NTS_vector(i) = mol_P3G_NTS;
        percent_NTS_vector_minus(i) = percent_NTS;
        
    else
        % If no solutions, assign NaN to indicate missing value
        mol_P3G_NTS_vector(i) = 0;
        percent_NTS_vector_minus(i) = 0;
    end
end

independentVariable = t_experimental;
dependentValues = percent_NTS_vector;
upperLimits = percent_NTS_vector_plus;
lowerLimits = percent_NTS_vector_minus;

% Calculate the error bars
errorBars = [dependentValues - lowerLimits, upperLimits - dependentValues];

% Define the value at which you want the horizontal line

% Define font sizes
labelFontSize = 12;   % Font size for axis labels
tickFontSize = 12;    % Font size for tick marks

% Define the RGB color triplet
scatterColor = [0.098, 0.612, 0.098]; % Corresponds to green color

% Scatter size
markerSize = 12;

% Create the main figure
figure;
set(gcf, 'Color', [1, 1, 1]); % Set figure background to white

% Main plot (figure 2 content)
mainAxes = axes;
scatter(mainAxes, independentVariable, dependentValues, markerSize, 'o', ...
    'MarkerFaceColor', scatterColor, 'MarkerEdgeColor', scatterColor); % Use RGB color
hold(mainAxes, 'on');
errorbar(mainAxes, independentVariable, dependentValues, errorBars(:,1), errorBars(:,2), ...
    'LineStyle', 'none', 'Color', 'k');

% Customize main plot
xlabel(mainAxes, 'time (min)', 'FontSize', labelFontSize);
ylabel(mainAxes, 'glycolytic flux through rubisco (%)', 'FontSize', labelFontSize);
%xlim(mainAxes, [min(independentVariable) - 1, max(independentVariable) + 1]);
xlim(mainAxes, [min(independentVariable) - 1, 140]);
ylim(mainAxes, [min(lowerLimits) - 1, max(upperLimits) + 1]);

% Customize main plot ticks
set(mainAxes, 'FontSize', tickFontSize);

hold(mainAxes, 'off');

% Create inset axes (figure 1 content)
insetAxes = axes('Position', [0.65, 0.3, 0.25, 0.25]); % [x, y, width, height] for inset position
scatter(insetAxes, independentVariable, dependentValues, markerSize, 'o', ...
    'MarkerFaceColor', scatterColor, 'MarkerEdgeColor', scatterColor); % Use RGB color
hold(insetAxes, 'on');
errorbar(insetAxes, independentVariable, dependentValues, errorBars(:,1), errorBars(:,2), ...
    'LineStyle', 'none', 'Color', 'k');

% Customize inset plot
xlim(insetAxes, [min(independentVariable) - 1, 25]);
ylim(insetAxes, [min(lowerLimits) - 1, 15]);
set(insetAxes, 'FontSize', tickFontSize);
xlabel(insetAxes, 'time (min)', 'FontSize', 10);
ylabel(insetAxes, 'glycolytic flux (%)', 'FontSize', 10);

% Add the black dashed horizontal line at y = y_horizontal_line
yline(insetAxes, y_horizontal_line, '--k'); % '--k' specifies a dashed black line

set(insetAxes, 'Color', [1, 1, 1]); % Set inset background to white

hold(insetAxes, 'off');

% Save the resulting figure as a TIFF file with a resolution of 200 dpi
print('flux_rubisco_labeling', '-dtiff', '-r200');

%% Read TIFF image
tiffImage = imread('input_metabolic_network.tif');
% Resize the image for better quality
scaleFactor = 0.5; % You can adjust this scale factor as needed
resizedImage = imresize(tiffImage, scaleFactor);
% Display the image with improved quality
figure('Name', 'TIFF Image', 'Position', [100, 100, size(resizedImage, 2), size(resizedImage, 1)]);
imshow(resizedImage);
%loading the .mat file containing the list of reactions and the coordinates
load input_reactions_and_positions
rxnID = findRxnIDs(model,reactions);

glucose_flux = FBAsolution.x(findRxnIDs(model,'GlycP'));

    for i=1:length(rxnID)    
        %fluxes(i)=FBAsolution.x(rxnID(i));
        fluxes(i) = (FBAsolution.x(rxnID(i))/glucose_flux);
        labels{i} = num2str(sprintf('%g',round(fluxes(i)*1000)/1000));
        text(axis_x(i),axis_y(i),labels(i),'FontSize',8,'Color','blue');
    end
    
% Save the resulting figure as a TIFF file with a resolution of 200 dpi
print('cartoon_labeling', '-dtiff', '-r200');
