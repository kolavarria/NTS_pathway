clc
clearvars
% Specify the file name to be deleted
fileToDelete = 'global_balance.txt';

% Check if the file exists before attempting to delete
if exist(fileToDelete, 'file')
        delete(fileToDelete);
else
        disp(['File "', fileToDelete, '" does not exist.']);
end

diary('global_balance.txt');



load Accumulibacter_anaerobic.mat

model = changeRxnBounds(model,{'TktA'},0,'b');
model = changeRxnBounds(model,{'ACS'},0,'l');
model = changeRxnBounds(model,{'GAPDH'},0,'l');
model = changeRxnBounds(model,{'GlycP'},1,'b');

FBAsolution = optimizeCbModel(model,'max');
disp ('--------------EMP-------------')
printFluxVector(model,FBAsolution.x,1,1);

ActiveRxnsFormula = printRxnFormula(model,model.rxns,false);
% Define the string to be replaced and the replacement string
stringToReplace = {'  '};
replacementString = {' '};

% Iterate through the array of strings and replace the target string
for i = 1:numel(ActiveRxnsFormula)
    ActiveRxnsFormula{i} = strrep(ActiveRxnsFormula{i}, stringToReplace, replacementString);
end

% Define the string to be replaced and the replacement string
stringToReplace = {'->'};
replacementString = {'<=>'};

% Iterate through the array of strings and replace the target string
for i = 1:numel(ActiveRxnsFormula)
    ActiveRxnsFormula{i} = strrep(ActiveRxnsFormula{i}, stringToReplace, replacementString);
end

T1 = table(model.rxns,ActiveRxnsFormula,FBAsolution.x,'VariableNames',{'Reaction Name' 'Reaction Formula' 'Relative Flux'});

writetable(T1, 'EMP_full.csv');

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

rxns_to_remove = [list_exchange,list_zero_flux,'ADPPPT','AMPPPT','PHBsyn'];
%rxns_to_remove = [list_exchange,list_zero_flux];

% Identify rows to remove
rows_to_remove = ismember(T1.("Reaction Name"), rxns_to_remove);

% Remove rows from the table
T1(rows_to_remove, :) = [];

order = {'GlycP';'GPM';'PGI';'PPiPFK';'FBPald';'TPI';'GAPDH';'PGK';'PGM';'ENO';'PYK';'PDH';'ACS';'PPiase';'Thio';'AAR'}; % the desired order    
[~, idx] = ismember(T1.("Reaction Name"), order);
[~, sortorder] = sort(idx);
T1 = T1(sortorder,:);
writetable(T1, 'EMP.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
load Accumulibacter_anaerobic.mat

model = changeRxnBounds(model,{'GAPDH'},0,'b');
model = changeRxnBounds(model,{'ACS'},0,'l');
model = changeRxnBounds(model,{'GlycP'},1,'b');

FBAsolution = optimizeCbModel(model,'max');
disp ('------------NTS-----------------')
printFluxVector(model,FBAsolution.x,1,1);

%diary off;

ActiveRxnsFormula = printRxnFormula(model,model.rxns,false);
% Define the string to be replaced and the replacement string
stringToReplace = {'  '};
replacementString = {' '};

% Iterate through the array of strings and replace the target string
for i = 1:numel(ActiveRxnsFormula)
    ActiveRxnsFormula{i} = strrep(ActiveRxnsFormula{i}, stringToReplace, replacementString);
end

% Define the string to be replaced and the replacement string
stringToReplace = {'->'};
replacementString = {'<=>'};

% Iterate through the array of strings and replace the target string
for i = 1:numel(ActiveRxnsFormula)
    ActiveRxnsFormula{i} = strrep(ActiveRxnsFormula{i}, stringToReplace, replacementString);
end

T1 = table(model.rxns,ActiveRxnsFormula,FBAsolution.x,'VariableNames',{'Reaction Name' 'Reaction Formula' 'Relative Flux'});
writetable(T1, 'NTS_full.csv');

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

rxns_to_remove = [list_exchange,list_zero_flux,'ADPPPT','AMPPPT','PHBsyn'];
%rxns_to_remove = [list_exchange,list_zero_flux];

% Identify rows to remove
rows_to_remove = ismember(T1.("Reaction Name"), rxns_to_remove);

% Remove rows from the table
T1(rows_to_remove, :) = [];

order = {'GlycP';'GPM';'PGI';'PPiPFK';'FBPald';'TPI';'TktA';'SBPald';'SBPase';'TktB';'RibE';'RibI';'RbuK';'RbuCO';'PGM';'ENO';'PYK';'PDH';'ACS';'PPiase';'Thio';'AAR'}; % the desired order    
[~, idx] = ismember(T1.("Reaction Name"), order);
[~, sortorder] = sort(idx);
T1 = T1(sortorder,:);
writetable(T1, 'NTS.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
load Accumulibacter_anaerobic.mat

model = changeRxnBounds(model,{'GAPDH'},0,'b');
model = changeRxnBounds(model,{'ACS'},0,'l');
model = changeRxnBounds(model,{'GlycP'},1,'b');
model = removeRxns(model,{'SBPald','SBPase'});
model = addReaction(model,'TALA','F6P + E4P <=> G3P + S7P');

FBAsolution = optimizeCbModel(model,'max');
disp ('------------RS-----------------')
printFluxVector(model,FBAsolution.x,1,1);

diary off;

ActiveRxnsFormula = printRxnFormula(model,model.rxns,false);
% Define the string to be replaced and the replacement string
stringToReplace = {'  '};
replacementString = {' '};

% Iterate through the array of strings and replace the target string
for i = 1:numel(ActiveRxnsFormula)
    ActiveRxnsFormula{i} = strrep(ActiveRxnsFormula{i}, stringToReplace, replacementString);
end

% Define the string to be replaced and the replacement string
stringToReplace = {'->'};
replacementString = {'<=>'};

% Iterate through the array of strings and replace the target string
for i = 1:numel(ActiveRxnsFormula)
    ActiveRxnsFormula{i} = strrep(ActiveRxnsFormula{i}, stringToReplace, replacementString);
end

T1 = table(model.rxns,ActiveRxnsFormula,FBAsolution.x,'VariableNames',{'Reaction Name' 'Reaction Formula' 'Relative Flux'});
writetable(T1, 'RS_full.csv');

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

rxns_to_remove = [list_exchange,list_zero_flux,'ADPPPT','AMPPPT','PHBsyn'];
%rxns_to_remove = [list_exchange,list_zero_flux];

% Identify rows to remove
rows_to_remove = ismember(T1.("Reaction Name"), rxns_to_remove);

% Remove rows from the table
T1(rows_to_remove, :) = [];

order = {'GlycP';'GPM';'PGI';'PPiPFK';'FBPald';'TPI';'TktA';'TALA';'TktB';'RibE';'RibI';'RbuK';'RbuCO';'PGM';'ENO';'PYK';'PDH';'ACS';'PPiase';'Thio';'AAR'}; % the desired order    
[~, idx] = ismember(T1.("Reaction Name"), order);
[~, sortorder] = sort(idx);
T1 = T1(sortorder,:);
writetable(T1, 'RS.csv');