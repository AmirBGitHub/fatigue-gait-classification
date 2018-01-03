function R = DCM(Data_in,fs)

% Feedback corrected DCM 

% this code uses direction cosine matrix and updates it to find the 
% current orientation using gravitational part of accelerometer signal
% and magnetometer data

acc_data = Data_in(:,1:3);
gyro_data = Data_in(:,4:6);
mag_data = Data_in(:,7:9);

duration = (length(Data_in)-1)/fs;
dt = 1/fs;
t = (0:dt:duration)';

% % 4rd order high-pass butterworth filter
% % with cutoff frequency 0.25hz
% fc = 0.25;
% [b,a] = butter(4,fc,'high');
% 
% % filter data
% acc_lin(:,1) = filter(b,a,acc_data(:,1));
% acc_lin(:,2) = filter(b,a,acc_data(:,2));
% acc_lin(:,3) = filter(b,a,acc_data(:,3));

%acc_g = acc_data + acc_lin;
acc_g = repmat([0 -9.81 0],length(acc_data),1);

R = zeros(3,3,length(t));
R(:,:,1) = [1 0 0; 0 1 0; 0 0 1];  % initial DCM
wICorrection = 0;
Wrp = 0.5;
Wy = 0.5;
KP = 0.2;
KI = 0.2;
wCorrection = [0 0 0];
wGyro = gyro_data*(pi/180);  % Gyro rate in rad/s

    for i = 1:length(t)-1

        w(i,:) = wGyro(i+1,:) + wCorrection;

        dTh_x = w(i,1)*dt;
        dTh_y = w(i,2)*dt;
        dTh_z = w(i,3)*dt;

        R(:,:,i+1) = R(:,:,i)*[1 -dTh_z dTh_y; dTh_z 1 -dTh_x; -dTh_y dTh_x 1];

        X = R(1,:,i+1);
        Y = R(2,:,i+1);

        err = dot(X,Y);

        X_orth = X - 0.5*err*Y;
        Y_orth = Y - 0.5*err*X;
        Z_orth = cross(X_orth,Y_orth);

        X_norm = 0.5*(3-dot(X_orth,X_orth))*X_orth;
        Y_norm = 0.5*(3-dot(Y_orth,Y_orth))*Y_orth;
        Z_norm = 0.5*(3-dot(Z_orth,Z_orth))*Z_orth;

        R(:,:,i+1) = [X_norm; Y_norm; Z_norm];

        YawCorrectionPlane = cross(R(1,:,i+1),mag_data(i+1,:));
        RollPitchCorrectionPlane = cross(R(3,:,i+1),acc_g(i+1,:));

        TotalCorrection = Wrp*RollPitchCorrectionPlane + Wy*YawCorrectionPlane;
        %TotalCorrection = Wrp*RollPitchCorrectionPlane;

        wPCorrection = KP*TotalCorrection;
        wICorrection = wICorrection + KI*dt*TotalCorrection;

        wCorrection = wPCorrection + wICorrection;
    end

end
