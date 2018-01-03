function [Seg_points, State, M_i_k, StepTime] = Process_and_Segment(Data_in)


    %% Read Data

    % This code imports the acceleration data from a spreadsheet and uses the
    % butterworth function to create a high pass filter to get rid of gravity.

    % accelerometer data reading
    acc_data = Data_in(:,1:3);
    gyro_data = Data_in(:,4:6);

    duration = (length(Data_in)-1)/51.2;
    dt = 1/51.2;
    t = (0:dt:duration)';
    
    %% Kalman filter 

    acc_g = repmat([0 -9.81 0],length(acc_data),1);

    [State, Rot] = IMU_Kalman(Data_in, acc_g); 
    filtered_angle = State(:,7:9);      % filtered angle
    unbiased_gyro = State(:,10:12);     % unbiased gyro
    
    % finds the quaternions of the rotation in order to transfering
    % body frame set of data to global frame
    
    acc_glob = Quaternion(acc_data,Rot);
    gyro_glob = Quaternion(gyro_data,Rot);
    

    acc_lin = acc_glob + acc_g; 
    Data_in_glob = [acc_glob gyro_glob];

    %% Segmentation and Velocity calculation

    
    % Acceleration 

    Acceleration_magn = sqrt(sum(acc_lin.^2,2));

    % Stance points detection
    
    % Acceleration filter parameters
    n=4;
    fc=4;
    Fs=35;
    Wn = (2/Fs)*fc;
    [b,a]=butter(n,Wn,'low');

    Acceleration_filt = filter(b,a,acc_lin);
    Acceleration_filt_magn = filter(b,a,sqrt(sum(Acceleration_filt.^2,2)));
    
    % Segmentation based on peaks in acc data
    Seg_points = [];
    Peak_points = [];
    n = 1;
%     while n<length(t)-55
%         if max(Acceleration_filt_magn(n+25:n+50))/max(Acceleration_filt_magn(n:n+24)) > 1.2 && max(Acceleration_filt_magn(n:n+50)) > 5    
%             [~,I1] = max(Acceleration_filt_magn(n+20:n+55));
%             Peak_points = [Peak_points; n+I1+19];
%             n = n + 50;
%         else
%             n = n + 1;
%         end
%     end
    
    while n<length(t)-30
        if max(Acceleration_filt_magn(n+20:n+30))/max(Acceleration_filt_magn(n:n+20)) > 1.2 && max(Acceleration_filt_magn(n:n+30)) > 5    
            [~,I1] = max(Acceleration_filt_magn(n:n+30));
            Peak_points = [Peak_points; n+I1];
            n = n + 25;
        else
            n = n + 1;
        end
    end
    
   
    
    % finding min points for segmentation
    for i = 1:length(Peak_points)-1
        [~,I2] = min(Acceleration_filt_magn(Peak_points(i):Peak_points(i)+10)); 
        Seg_points = [Seg_points; Peak_points(i)+I2];
    end
    
   % Finding stance phases
    Y_rest = [];
    k = 1;
    for i = 1:length(Peak_points)-1
        if t(Peak_points(i)) + 0.5 < t(end)
            Stance_point_s = Peak_points(i);
            while Acceleration_filt_magn(Stance_point_s) > 3
                Stance_point_s = Stance_point_s + 1;
            end
            Stance_point_e = Stance_point_s + 10;
            Y_rest = [Y_rest;(Stance_point_s:Stance_point_e)'];
            Y_rest_segs{k,1} = (Stance_point_s:Stance_point_e);
            k = k + 1;
        end
    end
    
    % Making segments windows
    k = 1;
    for j = 2:length(Seg_points)
        S_s(k) = Seg_points(j-1);
        S_e(k) = Seg_points(j)-1;
        k = k + 1;
    end     

    axis = [1 2 3];

    for direc = axis

        % Velocity calculation

        Vel{1,direc} = zeros(length(t),1);
        Res_acc = zeros(length(S_e),1);
        for i = 1:length(S_s)-1
            Res_acc(i,1) = trapz(t(S_s(i):S_e(i)),Acceleration_filt(S_s(i):S_e(i),direc))/range(t(S_s(i):S_e(i))); 
            for j = S_s(i):S_e(i)
                Vel{1,direc}(j,:) = trapz(t(S_s(i):j+1),Acceleration_filt(S_s(i):j+1,direc)) - Res_acc(i,1)*range(t(S_s(i):j+1));
            end
        end

        % Position calculation

        for i = 1:length(t)-1
            Pos{1,direc}(i+1,:) = trapz(t(1:i+1),(Vel{1,direc}(1:i+1)));
        end

        % Jerk calculation

        for i = 1:length(t)-1
            Jrk{1,direc}(i,:) = diff(Acceleration_filt(i:i+1,direc))/dt;
        end

        

%         Axes = {'X','Y','Z'};
% 
%         figure;
%         subplot(3,2,1);
%         plot(t,acc_lin(:,direc),'r')
%         title(['Unfiltered Acceleration ', Axes(direc)]);
%         xlabel('Time(s)');
%         ylabel('Acceleration (m/s^2)');
%         subplot(3,2,3);
%         plot(t,Acceleration_filt(:,direc),'r')
%         hold on
%         plot(t(locs),Acceleration_filt(locs,direc),'.g')
%         hold on
%         plot(t(Seg_ponits),Acceleration_filt(Seg_ponits,direc),'*b')
%         title(['Filtered Acceleration ', Axes(direc)]);
%         xlabel('Time(s)');
%         ylabel('Acceleration (m/s^2)');
% 
%         subplot(3,2,2);
%         plot(t,Vel{1,direc},'r')
%         hold on
%         plot(t(Seg_ponits),Vel{1,direc}(Seg_ponits),'*b')
%         title(['Velocity ', Axes(direc)]);
%         xlabel('Time(s)');
%         ylabel('Velocity (m/s)');
% 
%         subplot(3,2,4);
%         plot(t(2:end),Jrk{1,direc},'r')
%         hold on
%         plot(t(Seg_ponits),Jrk{1,direc}(Seg_ponits),'*b')
%         title(['Jerk ', Axes(direc)]);
%         xlabel('Time(s)');
%         ylabel('Jerk (m/s^3)');
% 
%         subplot(3,2,6);        
%         plot(t,Pos{1,direc},'r')
%         title(['Position ', Axes(direc)]);
%         xlabel('Time(s)');
%         ylabel('Position (m)');
%         

    end
    

    
%     figure;
%     subplot(2,1,1);        
%     plot(t,Acceleration_filt_magn,'r')
%     hold on
%     plot(t(Seg_points),Acceleration_filt_magn(Seg_points),'.g')
%     %plot(t(Y_rest),0,'.k')
%     title('Acceleration Magnitude');
%     xlabel('Time(s)');
%     ylabel('Acceleration (m/s^2)');
% 
%     Acceleration_filt_magn(Y_rest) = 0;
%     subplot(2,1,2);        
%     plot(t,Acceleration_filt_magn,'r')
%     hold on
%     plot(t(Seg_points),Acceleration_filt_magn(Seg_points),'.g')
%     plot(t(Y_rest),0,'.k')
%     title('Acceleration Magnitude');
%     xlabel('Time(s)');
%     ylabel('Acceleration (m/s^2)');
    
    
    % Metrics Segmentation
    
    Velocity = sqrt(Vel{1,1}.^2 + Vel{1,2}.^2 + Vel{1,3}.^2);
    Jerk = sqrt(Jrk{1,1}.^2 + Jrk{1,2}.^2 + Jrk{1,3}.^2);

    m = [Pos{1,1}(2:end), Pos{1,2}(2:end), t(2:end), Velocity(2:end), t(2:end), Acceleration_filt_magn(2:end), t(2:end), Jerk, filtered_angle(2:end,1), filtered_angle(2:end,2), filtered_angle(2:end,1), unbiased_gyro(2:end,1), filtered_angle(2:end,2), unbiased_gyro(2:end,2),filtered_angle(2:end,3), unbiased_gyro(2:end,3)];
    
    StepTime = diff(m(Seg_points,3));
    SegStep = diff(Seg_points);
    
    for i = 1:size(m,2)
        k = 1;
        for j = 1:length(Seg_points)-1
            if SegStep(j)<90 && SegStep(j)>4
%                 Making the data start from zero in each segment window
                M_i_k{k,i} = m(S_s(j):S_e(j),i) - m(S_s(j),i);         
%                 M_i_k{k,i} = m(S_s(j):S_e(j),i);
                k = k + 1;
            end
        end
    end
    
        % rotate and reflect the position data to get the actual values
    for i = 1:length(M_i_k(:,1))
        % define the x- and y-data for the original line we would like to rotate
        x = M_i_k{i,1}';
        y = M_i_k{i,2}';
        % create a matrix of these points, which will be useful in future calculations
        v = [x;y];
        % choose a point which will be the center of rotation
        x_center = x(1);
        y_center = y(1);
        % create a matrix which will be used later in calculations
        center = repmat([x_center; y_center], 1, length(x));
        % define a degree counter-clockwise rotation matrix
        theta = -atan2(M_i_k{i,2}(end)-M_i_k{i,2}(1), M_i_k{i,1}(end)-M_i_k{i,1}(1));       % rotation angle
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        % do the rotation...
        s = v - center;     % shift points in the plane so that the center of rotation is at the origin
        so = R*s;           % apply the rotation about the origin
        vo = so + center;   % shift again so the origin goes back to the desired center of rotation
        M_i_k{i,1} = vo(1,:)';
        M_i_k{i,2} = -vo(2,:)';
    end

%     % Makinge the same length segment windows    
%     maxLength = 50;
%     
%     for k = 1:length(M_i_k(:,1))
%         if length(M_i_k{k,1}) > maxLength
%             maxLength = length(M_i_k{k,1});
%         end 
%     end
%     for k = 1:length(M_i_k(:,1))
%         for i = 1:16
%             M_i_k{k,i}(length(M_i_k{k,i}):maxLength) = M_i_k{k,i}(end);
%         end
%     end
%     
%  
% %     figure;
% %     
% %     plot(m(Peak_points-1,2*CombinationSelection-1), m(Peak_points-1,2*CombinationSelection), 'ro')
% %     hold on
% %     plot(m(Seg_points-1,2*CombinationSelection-1), m(Seg_points-1,2*CombinationSelection), 'b*')
% %     hold on
% %     plot(m(:,2*CombinationSelection-1), m(:,2*CombinationSelection), 'Color', [1 0.4 0])
% %     hold on
% % 
% %     title(['Segmented profile of ', DataMode]);
% %     xlabel(MotionComponent(2*CombinationSelection-1));
% %     ylabel(MotionComponent(2*CombinationSelection));

end









  