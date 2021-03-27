clear
clc

% reading GS
[~,RXNs,~] = xlsread('CompetitiveInhibitors.xlsx','K1:K189');
[~,REGs,~] = xlsread('CompetitiveInhibitors.xlsx','I1:I189');

% reading fingerprints for molecules
[~,REGs_CHEBI,~] = xlsread('CompetitiveInhibitors.xlsx','Regulators','A1:A454');
REGs_CID = xlsread('CompetitiveInhibitors.xlsx','Regulators','B1:B454');
[~,REGs_Names,~] = xlsread('CompetitiveInhibitors.xlsx','Regulators','C1:C454');

FP_Subs = xlsread('CompetitiveInhibitors.xlsx','FP_Substrates','A1:A134');
FP_Regs = xlsread('CompetitiveInhibitors.xlsx','FP_Regulators','A1:A408');

FP_Subs_Vals = xlsread('CompetitiveInhibitors.xlsx','FP_Substrates','B1:AMK134');
FP_Regs_Vals = xlsread('CompetitiveInhibitors.xlsx','FP_Regulators','B1:AMK408');

emptyCells_rxn = cellfun(@isempty,RXNs);
emptyCells_reg = cellfun(@isempty,REGs);
RXNs(find(emptyCells_reg + emptyCells_rxn == 1)) = [];
REGs(find(emptyCells_reg + emptyCells_rxn == 1)) = [];

% the metabolic model
model = readCbModel;

RXNs_unique = unique(RXNs);
REGs_unique = unique(REGs);

RXNs_index_in_unique = zeros(length(RXNs),1);
for r = 1:length(RXNs)
    for r_u = 1:length(RXNs_unique)
        if strcmp(RXNs_unique{r_u},RXNs{r})
            RXNs_index_in_unique(r) = r_u;
            break
        end
    end
    
end

REGs_index_in_unique = zeros(length(REGs),1);
for r = 1:length(REGs)
    for r_u = 1:length(REGs_unique)
        if strcmp(REGs_unique{r_u},REGs{r})
            REGs_index_in_unique(r) = r_u;
            break
        end
    end
end

% Constructing the gold standard network
GS = zeros(length(RXNs_unique),length(REGs_unique));
for r = 1:length(RXNs_index_in_unique)
    GS(RXNs_index_in_unique(r),REGs_index_in_unique(r)) = 1;
end


% Positive features construction
inx_positive = 0;

[r_pos,c_pos] = find(GS == 1);

for r = 1:length(r_pos)
    rxn = RXNs_unique(r_pos(r));
    SUBs_in_model = [];
    for rr = 1:length(model.rxns)
        if strcmpi(rxn,model.rxns{rr})
            SUBs_in_model = find(model.S(:,rr) < 0);
            if model.lb(rr) < 0
                SUBs_in_model = union(SUBs_in_model,(find(model.S(:,rr) > 0)));
            end
            break
        end
    end
    
    SUBs_in_FP = [];
    ss = 0;
    for s = 1:length(SUBs_in_model)
        if ~isempty(find(FP_Subs == SUBs_in_model(s)))
            ss = ss + 1;
            SUBs_in_FP{ss} = find(FP_Subs == SUBs_in_model(s));
        end
    end
    
    if length(SUBs_in_FP) > 0
        Inx_r = strfind(REGs_CHEBI,REGs_unique{c_pos(r)});
        Index_r = find(not(cellfun('isempty',Inx_r)));
        if ~isempty(Index_r)
            if ~isempty(find(FP_Regs == REGs_CID(Index_r)))
                
                TM_temp = zeros(length(SUBs_in_FP),1);
                for i = 1:length(SUBs_in_FP)
                    FP_s = FP_Subs_Vals(SUBs_in_FP{i},:);
                    FP_r = FP_Regs_Vals(find(FP_Regs == REGs_CID(Index_r)),:);
                    
                    N_s = sum(FP_s);
                    N_r = sum(FP_r);
                    N_sr = sum(FP_s .* FP_r);
                    
                    
                    TM_temp(i) = 1 - (N_sr/(N_s + N_r - N_sr + eps));
                end
                [~,ii] = min(TM_temp);
                FP_s = FP_Subs_Vals(SUBs_in_FP{ii},:);
                FP_r = FP_Regs_Vals(find(FP_Regs == REGs_CID(Index_r)),:);
                
                inx_positive = inx_positive + 1;
                Features_positive(inx_positive,:) = [FP_s,FP_r];
                POSs(inx_positive,1) = r;

                
            end
        end
    end
end

% Negative features construction
inx_negative = 0;

[r_neg,c_neg] = find(GS == 0);

for r = 1:length(r_neg)
    rxn = RXNs_unique(r_neg(r));
    SUBs_in_model = [];
    for rr = 1:length(model.rxns)
        if strcmpi(rxn,model.rxns{rr})
            SUBs_in_model = find(model.S(:,rr) < 0);
            if model.lb(rr) < 0
                SUBs_in_model = union(SUBs_in_model,(find(model.S(:,rr) > 0)));
            end
            break
        end
    end
    
    SUBs_in_FP = [];
    ss = 0;
    for s = 1:length(SUBs_in_model)
        if ~isempty(find(FP_Subs == SUBs_in_model(s)))
            ss = ss + 1;
            SUBs_in_FP{ss} = find(FP_Subs == SUBs_in_model(s));
        end
    end
    
    if length(SUBs_in_FP) > 0
        Inx_r = strfind(REGs_CHEBI,REGs_unique{c_neg(r)});
        Index_r = find(not(cellfun('isempty',Inx_r)));
        if ~isempty(Index_r)
            if ~isempty(find(FP_Regs == REGs_CID(Index_r)))
                
                TM_temp = zeros(length(SUBs_in_FP),1);
                
                for i = 1:length(SUBs_in_FP)
                    FP_s = FP_Subs_Vals(SUBs_in_FP{i},:);
                    FP_r = FP_Regs_Vals(find(FP_Regs == REGs_CID(Index_r)),:);
                    
                    N_s = sum(FP_s);
                    N_r = sum(FP_r);
                    N_sr = sum(FP_s .* FP_r);
                    
                    
                    TM_temp(i) = 1 - (N_sr/(N_s + N_r - N_sr + eps));
                end
                [~,ii] = min(TM_temp);
                FP_s = FP_Subs_Vals(SUBs_in_FP{ii},:);
                FP_r = FP_Regs_Vals(find(FP_Regs == REGs_CID(Index_r)),:);
                
                inx_negative = inx_negative + 1;
                Features_negative(inx_negative,:) = [FP_s,FP_r];
                NEGs(inx_negative,1) = r;
                
            end
        end
    end
end


%% several SVMs for predicting negative labels
% The labels for negatives are already in the folder labels.mat, but here
% is where the procedure for choice of negatives.

% disp('Prediction of negative samples...     (takes a while)');
% 
% rand_num_Pos = randperm(size(Features_positive,1));
% rand_num_Neg = randperm(size(Features_negative,1));
% 
% STR = [num2str(floor(length(rand_num_Neg)/length(rand_num_Pos))),' SVMs are going to be trained for predicting negative samples.'];
% disp(STR);
% 
% for i = 1: floor(length(rand_num_Neg)/length(rand_num_Pos))
%     Start_train_neg(i) = 1 + length(rand_num_Pos) * (i - 1);
%     end_train_neg(i) = length(rand_num_Pos) * i;
% end
% 
% Labels = zeros(length(rand_num_Neg),1);
% 
% for i = 1:length(Start_train_neg)
%     
%     STR = ['SVM number ', num2str(i),' for predicting negative samples.'];
%     disp(STR);
%     
%     Data_Train = [Features_positive(rand_num_Pos(1:end),:) ; Features_negative(rand_num_Neg(Start_train_neg(i):end_train_neg(i)),:)];
%     
%     Group_Train = [ones(length(rand_num_Pos),1) ; zeros(length(rand_num_Pos),1)];
%     
%     SVMModel = fitcsvm(Data_Train,Group_Train,'Holdout',0.1,'Standardize',true,'KernelFunction','rbf',...
%         'KernelScale','auto');
%     CompactSVMModel = SVMModel.Trained{1}; % Extract trained, compact classifier
%     
%     range_test_neg = setdiff(1:length(rand_num_Neg),rand_num_Neg(Start_train_neg(i):end_train_neg(i)));
%     Data_Test_Neg = [Features_negative(range_test_neg,:)];
%     
%     [label,score] = predict(CompactSVMModel,Data_Test_Neg);
%     
%     Labels(range_test_neg) = label + Labels(range_test_neg);
% end

load Labels
Neg_labels = find(Labels == 0);


%% Final svm
Num_of_models = 10;
for s = 1:Num_of_models
    s
    RXN = [];
    Train_number = size(Features_positive,1);
   
    rand_Neg = randperm(size(Neg_labels,1));
    
    Data_Train = [Features_positive(1:Train_number,:) ; Features_negative(Neg_labels(rand_Neg(1:Train_number)),:)];
    
    Group_Train = [ones(Train_number,1) ; zeros(Train_number,1)];
    
    SVMModel_h = fitcsvm(Data_Train,Group_Train,'Standardize',true,'KernelFunction','rbf',...
        'KernelScale','auto','OptimizeHyperparameters','all', ...
        'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', ...
        'expected-improvement-plus','Holdout',0.1,'ShowPlots',false));
   
      % predicting all interactions
    inx_predict = 0;
    for r = 1:length(model.rxns)
        SUBs_in_model = [];
        SUBs_in_model = find(model.S(:,r) < 0);
        if model.lb(r) < 0
            SUBs_in_model = union(SUBs_in_model,(find(model.S(:,r) > 0)));
        end
        
        SUBs_in_FP = [];
        ss = 0;
        for su = 1:length(SUBs_in_model)
            if ~isempty(find(FP_Subs == SUBs_in_model(su)))
                ss = ss + 1;
                SUBs_in_FP{ss} = find(FP_Subs == SUBs_in_model(su));
            end
        end
        
        if length(SUBs_in_FP) > 0
            RXN = union(RXN,model.rxnNames(r)); 
            for reg = 1:length(REGs_unique)
                Inx_r = strfind(REGs_CHEBI,REGs_unique{reg});
                Index_r = find(not(cellfun('isempty',Inx_r)));
                if ~isempty(Index_r)
                    if ~isempty(find(FP_Regs == REGs_CID(Index_r)))
                        
                        TM_temp = zeros(length(SUBs_in_FP),1);
                        
                        for i = 1:length(SUBs_in_FP)
                            FP_s = FP_Subs_Vals(SUBs_in_FP{i},:);
                            FP_r = FP_Regs_Vals(find(FP_Regs == REGs_CID(Index_r)),:);
                            
                            N_s = sum(FP_s);
                            N_r = sum(FP_r);
                            N_sr = sum(FP_s .* FP_r);
                            
                            
                            TM_temp(i) = 1 - (N_sr/(N_s + N_r - N_sr + eps));
                        end
                        [~,ii] = min(TM_temp);
                        FP_s = FP_Subs_Vals(SUBs_in_FP{ii},:);
                        FP_r = FP_Regs_Vals(find(FP_Regs == REGs_CID(Index_r)),:);
                        
                        Feature_vector = [FP_s,FP_r];
                        
                        [label,~] = predict(SVMModel_h,Feature_vector);
                        inx_predict = inx_predict + 1;
                        Result{inx_predict,1} = model.rxnNames{r};
                        Result{inx_predict,2} = model.rxns{r};
                        Result{inx_predict,3} = model.grRules{r};
                        Result{inx_predict,4} = length(SUBs_in_FP);
                        Result{inx_predict,5} = REGs_Names{Index_r};
                        Result{inx_predict,6} = min(TM_temp); 
                        Result{inx_predict,6+s} = label;
                    end
                end
            end
        end
    end
    
end

Score = zeros(size(Result,1),1);
RESULT = zeros(size(Result,1),Num_of_models);

for i = 1:size(Result,1)
    for j = 1:Num_of_models
        RESULT(i,j) = Result{i,j+6};
    end
end

for i = 1:size(Result,1)
    Score(i) = sum(RESULT(i,:));
end

%% Predicted Netwok
Network = cell(size(Result,1),3);
Network(:,1) = Result(:,1);
Network(:,2) = Result(:,2);
Network(:,3) = num2cell(Score);

xlswrite('Predicte_Network.xlsx',Network);