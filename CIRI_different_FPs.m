clear
clc

% reading GS
[~,RXNs,~] = xlsread('CompetitiveInhibitors.xlsx','K1:K189');
[~,REGs,~] = xlsread('CompetitiveInhibitors.xlsx','I1:I189');

% reading fingerprints for molecules
[~,REGs_CHEBI,~] = xlsread('CompetitiveInhibitors.xlsx','Regulators','A1:A454');
REGs_CID = xlsread('CompetitiveInhibitors.xlsx','Regulators','B1:B454');

FP_Subs = xlsread('CompetitiveInhibitors.xlsx','FP_Substrates','A1:A134');
FP_Regs = xlsread('CompetitiveInhibitors.xlsx','FP_Regulators','A1:A408');

FP_Subs_Vals_1 = xlsread('CompetitiveInhibitors.xlsx','FP_Substrates','B1:AMK134');
FP_Regs_Vals_1 = xlsread('CompetitiveInhibitors.xlsx','FP_Regulators','B1:AMK408');


FP_Subs_Vals_2 = xlsread('CompetitiveInhibitors.xlsx','FP_Substrates_512','B1:SS134');
FP_Regs_Vals_2 = xlsread('CompetitiveInhibitors.xlsx','FP_Regulators_512','B1:SS408');


FP_Subs_Vals_3 = xlsread('CompetitiveInhibitors.xlsx','FP_Substrates_64','B1:BM134');
FP_Regs_Vals_3 = xlsread('CompetitiveInhibitors.xlsx','FP_Regulators_64','B1:BM408');

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
                    FP_s = FP_Subs_Vals_1(SUBs_in_FP{i},:);
                    FP_r = FP_Regs_Vals_1(find(FP_Regs == REGs_CID(Index_r)),:);
                    
                    N_s = sum(FP_s);
                    N_r = sum(FP_r);
                    N_sr = sum(FP_s .* FP_r);
                    
                    
                    TM_temp(i) = 1 - (N_sr/(N_s + N_r - N_sr + eps));
                end
                [~,ii] = min(TM_temp);
                FP_s = FP_Subs_Vals_1(SUBs_in_FP{ii},:);
                FP_r = FP_Regs_Vals_1(find(FP_Regs == REGs_CID(Index_r)),:);
                
                inx_positive = inx_positive + 1;
                Features_positive_1(inx_positive,:) = [FP_s,FP_r];
                
                TM_temp = zeros(length(SUBs_in_FP),1);
                for i = 1:length(SUBs_in_FP)
                    FP_s = FP_Subs_Vals_2(SUBs_in_FP{i},:);
                    FP_r = FP_Regs_Vals_2(find(FP_Regs == REGs_CID(Index_r)),:);
                    
                    N_s = sum(FP_s);
                    N_r = sum(FP_r);
                    N_sr = sum(FP_s .* FP_r);
                    
                    
                    TM_temp(i) = 1 - (N_sr/(N_s + N_r - N_sr + eps));
                end
                [~,ii] = min(TM_temp);
                FP_s = FP_Subs_Vals_2(SUBs_in_FP{ii},:);
                FP_r = FP_Regs_Vals_2(find(FP_Regs == REGs_CID(Index_r)),:);
                
                Features_positive_2(inx_positive,:) = [FP_s,FP_r];
                
                TM_temp = zeros(length(SUBs_in_FP),1);
                for i = 1:length(SUBs_in_FP)
                    FP_s = FP_Subs_Vals_3(SUBs_in_FP{i},:);
                    FP_r = FP_Regs_Vals_3(find(FP_Regs == REGs_CID(Index_r)),:);
                    
                    N_s = sum(FP_s);
                    N_r = sum(FP_r);
                    N_sr = sum(FP_s .* FP_r);
                    
                    
                    TM_temp(i) = 1 - (N_sr/(N_s + N_r - N_sr + eps));
                end
                [~,ii] = min(TM_temp);
                FP_s = FP_Subs_Vals_3(SUBs_in_FP{ii},:);
                FP_r = FP_Regs_Vals_3(find(FP_Regs == REGs_CID(Index_r)),:);
                
                Features_positive_3(inx_positive,:) = [FP_s,FP_r];
                
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
                for i = 1:length(SUBs_in_FP)
                    FP_s = FP_Subs_Vals_1(SUBs_in_FP{i},:);
                    FP_r = FP_Regs_Vals_1(find(FP_Regs == REGs_CID(Index_r)),:);
                    
                    inx_negative = inx_negative + 1;
                    
                    Features_negative_1(inx_negative,:) = [FP_s,FP_r];
                    
                    FP_s = FP_Subs_Vals_2(SUBs_in_FP{i},:);
                    FP_r = FP_Regs_Vals_2(find(FP_Regs == REGs_CID(Index_r)),:);
                    
                    Features_negative_2(inx_negative,:) = [FP_s,FP_r];
                    
                    FP_s = FP_Subs_Vals_3(SUBs_in_FP{i},:);
                    FP_r = FP_Regs_Vals_3(find(FP_Regs == REGs_CID(Index_r)),:);
                    
                    Features_negative_3(inx_negative,:) = [FP_s,FP_r];
                end
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
Train_number = floor(9 * size(Features_positive_1,1)/10);
Test_number = size(Features_positive_1,1) - Train_number;
num_of_training = 10;
accuracy_1 = zeros(num_of_training,1);
AUCsvm_1 = zeros(num_of_training,1);
PRsvm_1 = zeros(num_of_training,1);

accuracy_2 = zeros(num_of_training,1);
AUCsvm_2 = zeros(num_of_training,1);
PRsvm_2 = zeros(num_of_training,1);

accuracy_3 = zeros(num_of_training,1);
AUCsvm_3 = zeros(num_of_training,1);
PRsvm_3 = zeros(num_of_training,1);

i = 1;
while i <= num_of_training
    i
    rand_num_Pos = randperm(size(Features_positive_1,1));
    rand_Neg = randperm(size(Neg_labels,1));
    %  Training data is selected randomly:the training set is balanced- having the same number of positives
    %  and negatives-- This is done for all different fingerprints
    Data_Train_1 = [Features_positive_1(rand_num_Pos(1:Train_number),:) ; Features_negative_1(Neg_labels(rand_Neg(1:Train_number)),:)];
    Data_Train_2 = [Features_positive_2(rand_num_Pos(1:Train_number),:) ; Features_negative_2(Neg_labels(rand_Neg(1:Train_number)),:)];
    Data_Train_3 = [Features_positive_3(rand_num_Pos(1:Train_number),:) ; Features_negative_3(Neg_labels(rand_Neg(1:Train_number)),:)];
    
    Group_Train = [ones(Train_number,1) ; zeros(Train_number,1)];
    
    SVMModel_h_1 = fitcsvm(Data_Train_1,Group_Train,'Standardize',true,'KernelFunction','rbf',...
        'KernelScale','auto','OptimizeHyperparameters','all', ...
        'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', ...
        'expected-improvement-plus','Holdout',0.1,'ShowPlots',false));
    SVMModel_h_2 = fitcsvm(Data_Train_2,Group_Train,'Standardize',true,'KernelFunction','rbf',...
        'KernelScale','auto','OptimizeHyperparameters','all', ...
        'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', ...
        'expected-improvement-plus','Holdout',0.1,'ShowPlots',false));
    SVMModel_h_3 = fitcsvm(Data_Train_3,Group_Train,'Standardize',true,'KernelFunction','rbf',...
        'KernelScale','auto','OptimizeHyperparameters','all', ...
        'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', ...
        'expected-improvement-plus','Holdout',0.1,'ShowPlots',false));
    
    
    % Test
    Data_Test_1 = [Features_positive_1(rand_num_Pos(Train_number+1:end),:) ; Features_negative_1(Neg_labels(rand_Neg(Train_number+1:size(Features_positive_1,1))),:)];
    Data_Test_2 = [Features_positive_2(rand_num_Pos(Train_number+1:end),:) ; Features_negative_2(Neg_labels(rand_Neg(Train_number+1:size(Features_positive_2,1))),:)];
    Data_Test_3 = [Features_positive_3(rand_num_Pos(Train_number+1:end),:) ; Features_negative_3(Neg_labels(rand_Neg(Train_number+1:size(Features_positive_3,1))),:)];
    Test_number = size(Features_positive_1,1) - Train_number;
    Group_Test = [ones(Test_number,1) ; zeros(Test_number,1)];
    
    
    [label_1,score_1] = predict(SVMModel_h_1,Data_Test_1);
    [label_2,score_2] = predict(SVMModel_h_2,Data_Test_2);
    [label_3,score_3] = predict(SVMModel_h_3,Data_Test_3);
    
    accuracy_1(i,1) = sum(predict(SVMModel_h_1, Data_Test_1) == Group_Test)/length(Group_Test)*100;
    accuracy_2(i,1) = sum(predict(SVMModel_h_2, Data_Test_2) == Group_Test)/length(Group_Test)*100;
    accuracy_3(i,1) = sum(predict(SVMModel_h_3, Data_Test_3) == Group_Test)/length(Group_Test)*100;
    
    [X_AUC_1,Y_AUC_1,Tsvm,AUCsvm_1(i,1)] = perfcurve(logical(Group_Test),score_1(:,logical(SVMModel_h_1.ClassNames)),'true');
    [X_PR_1,Y_PR_1,~,PRsvm_1(i,1)] = perfcurve(logical(Group_Test),score_1(:,logical(SVMModel_h_1.ClassNames)),'true','xCrit', 'reca', 'yCrit', 'prec');
    
    [X_AUC_2,Y_AUC_2,Tsvm,AUCsvm_2(i,1)] = perfcurve(logical(Group_Test),score_2(:,logical(SVMModel_h_2.ClassNames)),'true');
    [X_PR_2,Y_PR_2,~,PRsvm_2(i,1)] = perfcurve(logical(Group_Test),score_2(:,logical(SVMModel_h_2.ClassNames)),'true','xCrit', 'reca', 'yCrit', 'prec');
    
    [X_AUC_3,Y_AUC_3,Tsvm,AUCsvm_3(i,1)] = perfcurve(logical(Group_Test),score_3(:,logical(SVMModel_h_3.ClassNames)),'true');
    [X_PR_3,Y_PR_3,~,PRsvm_3(i,1)] = perfcurve(logical(Group_Test),score_3(:,logical(SVMModel_h_3.ClassNames)),'true','xCrit', 'reca', 'yCrit', 'prec');
    
    
    % This is just to see the cases with close AUC values
    if abs(AUCsvm_1(i,1) - AUCsvm_2(i,1)) <= 0.05
        if abs(AUCsvm_3(i,1) - AUCsvm_2(i,1)) <= 0.05
            break
        end
    end
    i = i + 1;

    
end
%% Display results for the classifier with difference between the FPs less than 0.05

figure(1)
hold on
t = linspace(realmin ( 'single' ),1);
plot(t,t,'--','Color',[0,0.45,0.74]);
plot(X_AUC_1,Y_AUC_1,'LineWidth',1.5,'Color',[0 0.4470 0.7410]);
plot(X_AUC_2,Y_AUC_2,'LineWidth',1.5,'Color',[0.4660 0.6740 0.1880]);
plot(X_AUC_3,Y_AUC_3,'LineWidth',1.5,'Color',[0.6350 0.0780 0.1840]	);

xlabel('False positive rate')
ylabel('True positive rate')
title('ROC')

figure(2)
hold on
t = linspace(realmin ( 'single' ),1);
plot(t,1-t,'--','Color',[0,0.45,0.74]);
plot(X_PR_1,Y_PR_1,'LineWidth',1.5,'Color',[0 0.4470 0.7410]);
plot(X_PR_2,Y_PR_2,'LineWidth',1.5,'Color',[0.4660 0.6740 0.1880]);
plot(X_PR_3,Y_PR_3,'LineWidth',1.5,'Color',[0.6350 0.0780 0.1840]);

xlabel('Recall')
ylabel('Precision')
title('PR')
