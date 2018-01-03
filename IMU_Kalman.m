function [State, R] = IMU_Kalman(imu_data, acc_g)

duration = (length(imu_data)-1)/51.2;
dt = 1/51.2;
t = (0:dt:duration)';

acc_data = imu_data(:,1:3);
gyro_data = imu_data(:,4:6);
gyro = gyro_data;

v_gyro = sqrt(0.3)*randn(length(t),3); % measurement noise for gyro,  variance = 0.3
v_acc = sqrt(0.3)*randn(length(t),3); % measurement noise for accellerometer, variance = 0.3

% Compute the angles computed by using only accelerometers of gyroscope
x = [atan(acc_data(:,1)./sqrt(sum([acc_data(:,2).^2 acc_data(:,3).^2],2)))*(180/pi), atan(acc_data(:,2)./sqrt(sum([acc_data(:,1).^2 acc_data(:,3).^2],2)))*(180/pi), atan(acc_data(:,3)./sqrt(sum([acc_data(:,1).^2 acc_data(:,2).^2],2)))*(180/pi)];
x_acc = x + v_acc; % Angle computed from accelerometer measurement
% x_acc = x;
x_gyro = cumsum(gyro*dt,1); % Angle computed by integration of gyro measurement


P = [1 0; 0 1];
R_angle = 0.3;

Q_angle = 0.05;
Q_gyro = 0.5;
Q = [Q_angle 0; 0 Q_gyro];

A = [0 -1; 0 0];
q_bias = [0 0 0]; % Initialize gyro bias
angle = [0 90 0]; % Initialize gyro angle
q_m = 0;
X = [0 90 0; 0 0 0];

    for i=1:length(t)

         % Gyro update 

         q_m = gyro(i,:);

         q = q_m - q_bias; % gyro bias removal

         Pdot = A*P + P*A' + Q;

         rate = q;

         angle = angle + q*dt;

         P = P + Pdot*dt;

         % Kalman (Accelerometer) update 

         C = [1 0];
         angle_err = x_acc(i,:)-angle;
         E = C*P*C' + R_angle;

         K = P*C'*inv(E);

         P = P - K*C*P;
         X = X + K * angle_err;
         x1(i,:) = X(1,:);
         R(:,:,i) = ang2orth(x1);
         x2(i,:) = X(2,:);
         angle = x1(i,:);
         q_bias = x2(i,:);

         x3(i,:) = q;  % unbiased gyro rate
    end

State = [x_acc x_gyro x1 x3];

x_acc = State(:,1:3);   % angle calculated from accelerometer
x_gyro = State(:,4:6);  % angle calculated from gyro 
filtered_angle = State(:,7:9);      % filtered angle
unbiased_gyro = State(:,10:12);     % unbiased gyro
rotation_mat = R;


% %Plot the state before using Kalman filter
% figure;
% subplot(3,4,1);
% plot(t,x_acc(:,1));
% xlabel('time(s)');
% ylabel('Theta x(t)');
% legend('acc angle');
% 
% subplot(3,4,5);
% plot(t,x_acc(:,2));
% xlabel('time(s)');
% ylabel('Theta y(t)');
% legend('acc angle');
% 
% subplot(3,4,9);
% plot(t,x_acc(:,3));
% xlabel('time(s)');
% ylabel('Theta z(t)');
% legend('acc angle');
% 
% subplot(3,4,2);
% plot(t,x_gyro(:,1));
% xlabel('time(s)');
% ylabel('Theta x(t)');
% legend('gyro angle');
% 
% subplot(3,4,6);
% plot(t,x_gyro(:,2));
% xlabel('time(s)');
% ylabel('Theta y(t)');
% legend('gyro angle');
% 
% subplot(3,4,10);
% plot(t,x_gyro(:,3));
% xlabel('time(s)');
% ylabel('Theta z(t)');
% legend('gyro angle');
% 
% % Plot the result using kalman filter
% subplot(3,4,3);
% plot(t,unbiased_gyro(:,1));
% xlabel('time(s)');
% ylabel('ThetaDot x(t)');
% legend('gyro rate unbiased');
% 
% subplot(3,4,7);
% plot(t,unbiased_gyro(:,2));
% xlabel('time(s)');
% ylabel('ThetaDot y(t)');
% legend('gyro rate unbiased');
% 
% subplot(3,4,11);
% plot(t,unbiased_gyro(:,3));
% xlabel('time(s)');
% ylabel('ThetaDot z(t)');
% legend('gyro rate unbiased');
% 
% subplot(3,4,4);
% plot(t,filtered_angle(:,1));
% xlabel('time(s)');
% ylabel('Theta x(t)');
% legend('kalman angle');
% 
% subplot(3,4,8);
% plot(t,filtered_angle(:,2));
% xlabel('time(s)');
% ylabel('Theta y(t)');
% legend('kalman angle');
% 
% subplot(3,4,12);
% plot(t,filtered_angle(:,3));
% xlabel('time(s)');
% ylabel('Theta z(t)');
% legend('kalman angle');

    function orthm = ang2orth(ang) 
        sa = sind(ang(2)); ca = cosd(ang(2)); 
        sb = sind(ang(1)); cb = cosd(ang(1)); 
        sc = sind(ang(3)); cc = cosd(ang(3)); 

        ra = [  ca,  sa,  0; ... 
               -sa,  ca,  0; ... 
                 0,   0,  1]; 
        rb = [  cb,  0,  sb; ... 
                 0,  1,  0; ... 
               -sb,  0,  cb]; 
        rc = [  1,   0,   0; ... 
                0,   cc, sc;... 
                0,  -sc, cc]; 
        orthm = rc*rb*ra;
    end

end