clear all
close all
clc

% Motion Component Selection

Feature = {'Position x, Position y';                % #1
           't, Velocity';                           % #2
           't, Acceleration';                       % #3
           't, Jerk';                               % #4
           'Theta x, Theta y';                      % #5
           'Theta x, Theta_dot x';                  % #6
           'Theta y, Theta_dot y';                  % #7
           'Theta z, Theta_dot z'};                 % #8


MotionComponent = {'Position x (m)', 'Position y (m)', 'time(s)', 'Velocity Magnitude (m/s)', 'time(s)', 'Acceleration Magnitude (m/s)', 'time(s)', 'Jerk Magnitude (m/s)', 'Filtered Angle x (degree)', 'Filtered Angle y (degree)', 'Filtered Angle x (degree)', 'Unbiased Gyro x (degree/s)', 'Filtered Angle y (degree)', 'Unbiased Gyro y (degree/s)', 'Filtered Angle z (degree)', 'Unbiased Gyro z (degree/s)'};

% load ScoreStepTimeCrossValSubj20.mat
% load TrueLabelCrossValSubj20.mat

for sbjct = 1:20
    
    %% Data Read, Process and Save 
    
    Participant = Subject(sbjct)
    m = sbjct;
     
    clear imu_data
    % Read Data 

    imu_read = importdata([pwd, '\Testdata_(', char(sbjct), ').xlsx']);
    imu_data = imu_read.data.Ankle;
%     imu_data = imu_read.data;

    Data_in_S = imu_data(:,2:7);
    Data_in_E = imu_data(:,13:18);
    Data_in_tS = imu_data(:,24:29);
    Data_in_tE = imu_data(:,35:40);

    % Data Processing and Segmentation
    CombinationSelection = 1;       
    [Seg_points_S, State_S, M_i_k_S, StepTime_S] = Process_and_Segment(Data_in_S);
    save(['MIK_S', char(Subject(sbjct))],'M_i_k_S')
    [Seg_points_E, State_E, M_i_k_E, StepTime_E] = Process_and_Segment(Data_in_E);
    save(['MIK_E', char(Subject(sbjct))],'M_i_k_E')
    [Seg_points_tS, State_tS, M_i_k_tS, StepTime_tS] = Process_and_Segment(Data_in_tS);
    save(['MIK_tS', char(Subject(sbjct))],'M_i_k_tS')
    [Seg_points_tE, State_tE, M_i_k_tE, StepTime_tE] = Process_and_Segment(Data_in_tE);
    save(['MIK_tE', char(Subject(sbjct))],'M_i_k_tE')

    %% Step Time Calculation
    
    MIK_S_load = load(['MIK_S', char(Subject(sbjct))]); 
    M_i_k_S = MIK_S_load.M_i_k_S(end-25:end-1,:); 
    for i = 1:length(M_i_k_S(:,1))
      StepTime_S(i,1) = range(M_i_k_S{i,3});
    end
    
    MIK_E_load = load(['MIK_E', char(Subject(sbjct))]); 
    M_i_k_E = MIK_E_load.M_i_k_E(end-25:end-1,:); 
    for i = 1:length(M_i_k_E(:,1))
      StepTime_E(i,1) = range(M_i_k_E{i,3});
    end
    
    MIK_tS_load = load(['MIK_tS', char(Subject(sbjct))]); 
    M_i_k_tS = MIK_tS_load.M_i_k_tS(end-25:end-1,:);
    
    MIK_tE_load = load(['MIK_tE', char(Subject(sbjct))]); 
    M_i_k_tE = MIK_tE_load.M_i_k_tE(end-25:end-1,:);
     

    %% Testing and Training Set Initialization

    Non_Fatigued_SamplCnt = size(M_i_k_S,1); 
    Fatigued_SamplCnt = size(M_i_k_E,1); 
    Non_Fatigued_TstCnt = size(M_i_k_tS,1);
    Fatigued_TstCnt = size(M_i_k_tE,1);

    Threshold_Class = Non_Fatigued_SamplCnt; 
    Threshold_Test = Non_Fatigued_TstCnt;

    TestingData = [M_i_k_tS; M_i_k_tE];   % first 40sec and last 40sec for training in 1$

    for i = 1:length(TestingData(:,1))
        StepTime_t(i,1) = range(TestingData{i,3});
    end

    TrainingData = [M_i_k_S; M_i_k_E];    % first 40sec and last 40sec for training in 1$

    % Generating random test sample 
    TestIndex = randperm(length(TestingData),50);
    
               
    for CombinationSelection = 1:8

        n = CombinationSelection;
        
        lNF = []; 
        lF = [];
        
        for trial = 1:50

            % Test data picking 
            l = trial;

            TestData = [TestingData{TestIndex(l),2*n-1} TestingData{TestIndex(l),2*n}];

            
            %% 1$ Recognizer

            % Distance based scoring in 1$ recognizer 
            
            [Index1, Score1] = OneDollar_Recognizer_Euclid(TestData,M_i_k_S,n,25); 
            [Index2, Score2] = OneDollar_Recognizer_Euclid(TestData,M_i_k_E,n,25);
            
            if TestIndex(l) <= Threshold_Test
                lNF = [lNF;l]; 
                TrueLabelAll(l,1,m) = -1;
                %ScoreStepTime(l,1:2,n) = [mean(Score1), abs(mean(StepTime_S(1:25))-StepTime_t(TestIndex(l)))];
                ScoreStepTime(l,1:2,n,m) = [mean(Score1), StepTime_t(TestIndex(l))];
            elseif TestIndex(l) > Threshold_Test
                lF = [lF;l]; 
                TrueLabelAll(l,1,m) = 1;
                ScoreStepTime(l,1:2,n,m) = [mean(Score2), StepTime_t(TestIndex(l))];
            end

                          
        end
        
        
%         %% Graphical Demonstration
% 
%         figure('Color',[1 1 1]); 
%         plot(ScoreStepTime(find(TrueLabelAll(:,1,1)==-1),2,n,m), ScoreStepTime(find(TrueLabelAll(:,1,m)==-1),1,n,m), 'k*', ScoreStepTime(find(TrueLabelAll(:,1,1)==1),2,n,m), ScoreStepTime(find(TrueLabelAll(:,1,m)==1),1,n,m), 'k+', 'linewidth',1.1, 'MarkerSize', 10)
%         hold on
%         hLeg = legend('Non-Fatigued','Fatigued');    
%         set(hLeg,'FontSize',14, 'interpreter', 'latex');
%         xlabel('Step Time (sec)','interpreter','latex','FontSize', 14,'FontWeight','bold');
%         ylabel('Euclidean Distance Based Score','interpreter','latex', 'FontSize',14,'FontWeight','bold');  


        %% Train and Cross Validate and Label Test Sample Observations of
        % SVM Classifiers
        
        for f = 1:5

            ftrain1 = ScoreStepTime([1:10*(f-1), 10*f+1:50],1,n,m);
            ftrain2 = ScoreStepTime([1:10*(f-1), 10*f+1:50],2,n,m);
            grouptrain = TrueLabelAll([1:10*(f-1), 10*f+1:50],1,m);
            ftest1 = ScoreStepTime(10*(f-1)+1:10*f,1,n,m); 
            ftest2 = ScoreStepTime(10*(f-1)+1:10*f,2,n,m); 
            grouptest = TrueLabelAll(10*(f-1)+1:10*f,1,m);


            % Train an SVM classifier. 
            CVSVMModel = fitcsvm([ftrain1 ftrain2],grouptrain,'KernelFunction','RBF','Standardize',true);

            XTest = [ftest1 ftest2]; 
            YTest = grouptest;
            
            [label,score] = predict(CVSVMModel,XTest);
            

%             table(YTest(1:10),label(1:10),'VariableNames',{'TrueFatigue','PredictedFatigue'})
            TrueLabel(:,1,f) = grouptest; 
            PredictedLabel(:,n,f) = label;
                                        
        end
                  
    end
    
    %% Classification Results
    
    for f = 1:5
        ClassifyLabels = [TrueLabel(:,1,f) NaN(10,1) PredictedLabel(:,:,f) NaN(10,1) 2*(sum(PredictedLabel(:,:,f),2)>0)-1]; 
        j = 1; 
        for i = [3 4 5 6 7 8 9 10 12]
            [C,order] = confusionmat(ClassifyLabels(:,1),ClassifyLabels(:,i));
            Accuracy(f,j) = 100*(C(1,1)+C(2,2))/(C(1,1)+C(2,2)+C(2,1)+C(1,2));
            Sensitivity(f,j) = 100*C(2,2)/(C(2,2)+C(2,1));
            Specificity(f,j) = 100*C(1,1)/(C(1,1)+C(1,2)); 
            j = j + 1;
        end
    end
                       
    AccuracyFolds(m,:) = mean(Accuracy); 
    SensitivityFolds(m,:) = mean(Sensitivity); 
    SpecificityFolds(m,:) = mean(Specificity);
   
end

AccuracyTotal = mean(AccuracyFolds); 
SensitivityTotal = mean(SensitivityFolds); 
SpecificityTotal = mean(SpecificityFolds);

    
    

