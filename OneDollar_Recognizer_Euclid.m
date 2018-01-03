function [index, score] = OneDollar_Recognizer_Euclid(TestData,TrainingData,MotionComponent,TrainSize)
    
%     M = TemplatesCtn; % number of templates
    m = MotionComponent; % the combination of motion components selected
    M = TrainSize;
    N = 64;
    Templates = zeros(N,2,M);
    
    %% Step 1 this part is removed beacuse we have made the same length data before
    
    for c = 1:M
%         Templates(:,:,c) = Resample(TrainingData{1,c},N);
        Templates(:,:,c) = Resample([TrainingData{c,2*m-1} TrainingData{c,2*m}],N);
    end
    Data = Resample(TestData,N);

    
    %% Step 2
    
%     for c = 1:M
%         Templates(:,:,c) = Rotate_To_Zero([TrainingData{c,2*m-1} TrainingData{c,2*m}]);
%     end
%     Data = Rotate_To_Zero(TestData);
    
    for c = 1:M
        Templates(:,:,c) = Rotate_To_Zero(Templates(:,:,c));
    end
    Data = Rotate_To_Zero(Data);
     
    %% Step 3
    
    Size = sqrt((max(Data(:,1)) - min(Data(:,1)))^2 + (max(Data(:,2)) - min(Data(:,2)))^2);
    for c = 1:M
        Templates(:,:,c) = Scale_To_Square(Templates(:,:,c),Size);
        Templates(:,:,c) = Translate_To_Origin(Templates(:,:,c));
    end
    Data = Scale_To_Square(Data,Size);
    Data = Translate_To_Origin(Data);
    
    %% Step 4
    
    [score,index] = Recognize(Data,Templates,Size,M);
    
    %% Functions
    
    function newPoints = Resample(Points,n)
        I = Path_Length(Points)/(n-1);
        D = 0;
        newPoints = Points(1,:);
        i = 2;
        while(i <= length(Points))
            d = pdist2(Points(i-1,:),Points(i,:));
            if (D+d)>=I
                q = Points(i-1,:) + ((I-D)/d)*(Points(i,:)-Points(i-1,:));
                newPoints = [newPoints ; q];
                Points = [Points(1:i-1,:) ; q ; Points(i:end,:)];
                D = 0;
            else
                D = D+d;
            end
            i = i+1;
        end
        if length(newPoints)~=n
            newPoints = [newPoints ; Points(end,:)];
        end
    end

    function d = Path_Length(A)
        d = 0;
        for i = 1:length(A)-1
            d = d + pdist2(A(i,:),A(i+1,:));
        end
    end

    function newPoints = Rotate_To_Zero(Points)
        C = mean(Points);
        Theta = atan2(C(2)-Points(1,2) , C(1)-Points(1,1));
        newPoints = Rotate_By(Points,-Theta);
    end

    function newPoints = Rotate_By(Points,Theta)
        C = mean(Points);
        newPoints = Points*0;
        for i=1:length(Points)
            qx = (Points(i,1)-C(1))*cos(Theta) - (Points(i,2)-C(2))*sin(Theta) + C(1);
            qy = (Points(i,1)-C(1))*sin(Theta) + (Points(i,2)-C(2))*cos(Theta) + C(2);
            newPoints(i,:) = [qx qy];
        end
    end

    function newPoints = Scale_To_Square(Points,Size)
        B_width =  max(Points(:,1)) - min(Points(:,1));
        B_height =  max(Points(:,2)) - min(Points(:,2));
        B = [B_width, B_height];  
        newPoints = Points*diag(Size./B);
    end

    function newPoints = Translate_To_Origin(Points)
       C = mean(Points);
       newPoints = Points - repmat(C,length(Points),1);
    end
        
    function [score, index] = Recognize(Points,templates,Size,M)
        b = inf;
        index = 0;
        Theta = pi/4;
        Theta_Delta = 2*pi/180;
        for i = 1:M
            T = templates(:,:,i);
            d = Distance_at_Best_Angle(Points,T,-Theta,Theta,Theta_Delta);
            %if d < b
                b = d;
                T_p = T;
                %index = i;
            %end
            score(i,1) = 1 - b/(0.5*sqrt(Size^2+Size^2));
        end
        [~,index] = max(score);
        %score = 1 - b/(0.5*sqrt(Size^2+Size^2));
    end

    function f = Distance_at_Best_Angle(Points,T,Theta_a,Theta_b,td)
        Phi = 0.5*(-1 + sqrt(5));
        x1 = Phi*Theta_a + (1-Phi)*Theta_b;
        f1 = Distance_at_Angle(Points,T,x1);
        x2 = (1-Phi)*Theta_a + Phi*Theta_b;
        f2 = Distance_at_Angle(Points,T,x2);
        while(abs(Theta_b-Theta_a)>td)
            if f1<f2
                Theta_b = x2;
                x2 = x1;
                f2 = f1;
                x1 = Phi*Theta_a + (1 - Phi)*Theta_b;
                f1 = Distance_at_Angle(Points,T,x1);
            else
                Theta_a = x1;
                x1 = x2;
                f1 = f2;
                x2 = (1 - Phi)*Theta_a + Phi*Theta_b;
                f2 = Distance_at_Angle(Points,T,x2);
            end
        end
        f = min(f1,f2);
    end

    function d = Distance_at_Angle(Points,T,Th)
        newPoints = Rotate_By(Points,Th);
        d = Path_Distance(newPoints,T);
    end

    function d_p = Path_Distance(A,B)
        d = 0;
        for i = 1:length(A)
            d = d + pdist2(A(i,:),B(i,:));
        end
        d_p = d/length(A);
    end
    
end
















