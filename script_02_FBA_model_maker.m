%Script to construct a .mat model
clc
clearvars
model = createModel();

model=addReaction(model,'GlycP','Mal4 + Pi -> G1P + Mal3');
model=addReaction(model,'GPM','G1P <=> G6P');
model=addReaction(model,'PGI','G6P <=> F6P');
model=addReaction(model,'PPiPFK','F6P + PPi -> FBP + Pi');
model=addReaction(model,'FBPald','FBP <=> DHAP + G3P');
model=addReaction(model,'TPI','DHAP <=> G3P');
model=addReaction(model,'GAPDH','G3P + NAD + Pi <=> BPG + NADH');
model=addReaction(model,'PGK','BPG + ADP <=> P3G + ATP');
model=addReaction(model,'PGM','P3G <=> P2G');
model=addReaction(model,'ENO','P2G <=> PEP + H2O');
model=addReaction(model,'PYK','PEP + ADP -> ATP + PYR');
model=addReaction(model,'PDH','PYR + NAD + CoA -> AcCoA + NADH + CO2');
model=addReaction(model,'ACS','Ace + CoA + ATP -> AcCoA + PPi + AMP');
model=addReaction(model,'Thio','AcCoA + AcCoA -> AACoA + CoA');
model=addReaction(model,'AAR','AACoA + NADH -> HBCoA + NAD');
model=addReaction(model,'PHBsyn','HBCoA + HB -> PHB + CoA');
model=addReaction(model,'TktA','F6P + G3P <=> X5P + E4P');
model=addReaction(model,'SBPald','DHAP + E4P <=> SBP');
model=addReaction(model,'SBPase','SBP + H2O -> S7P + Pi');
model=addReaction(model,'TktB','G3P + S7P <=> Ri5P + X5P');
model=addReaction(model,'RibE','X5P <=> Ru5P');
model=addReaction(model,'RibI','Ri5P <=> Ru5P');
model=addReaction(model,'RbuK','Ru5P + ATP -> RBP + ADP');
model=addReaction(model,'RbuCO','RBP + CO2 + H2O -> P3G + P3G');
model=addReaction(model,'ADPPPT','PolyPP + ADP -> ATP + PolyP');
model=addReaction(model,'AMPPPT','PolyPP + AMP -> ADP + PolyP');
model=addReaction(model,'PPiase','PPi + H2O -> Pi + Pi');

model = addExchangeRxn(model,{'PolyPP'},-50,50);
model = addExchangeRxn(model,{'PolyP'},-50,50);
model = addExchangeRxn(model,{'Mal4'},-5,50);
model = addExchangeRxn(model,{'Mal3'},-50,50);
model = addExchangeRxn(model,{'Ace'},-50,50);
model = addExchangeRxn(model,{'Pi'},0,50);
model = addExchangeRxn(model,{'CO2'},-50,50);
model = addExchangeRxn(model,{'HB'},-50,50);
model = addExchangeRxn(model,{'PHB'},-50,50);
model = addExchangeRxn(model,{'H2O'},-50,50);

model = changeObjective(model,'EX_PHB');

save Accumulibacter_anaerobic
clc

% disp ('----in silico model-----------')

% model = changeRxnBounds(model,{'GAPDH'},0,'b');
% model = changeRxnBounds(model,{'GlycP'},1,'b');

% model = removeRxns(model,{'SBPald','SBPase'});
% model = addReaction(model,'TALA','F6P + E4P <=> G3P + S7P');

RxnsFormula = printRxnFormula(model,model.rxns,false);
T0 = table(model.rxns,RxnsFormula,'VariableNames',{'Reaction Name' 'Reaction Formula'});
writetable(T0,'in_silico_model.csv');

% FBAsolution = optimizeCbModel(model,'max');
% disp ('------- q-rates -------------')
% printFluxVector(model, FBAsolution.x,1,1);

% disp ('--------- flux distributions --------------------')
% T1=table(model.rxns,RxnsFormula,FBAsolution.x,'VariableNames',{'Reaction Name' 'Reaction Formula' 'Fluxes'});
% disp(T1)
