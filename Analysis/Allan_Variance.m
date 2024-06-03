data = rosbag("C:/Users/visar/Downloads/LocationA.bag");
% rosbag info LocationA.bag;
readings = select(data,'Topic','/vectornav');
topics=readings.AvailableTopics;

msg = readMessages(readings,'DataFormat','struct');
%msg{1};
%writetable(struct2table(msg{1}), 'LocationB.csv');

for i = 1:length(msg)
    val = split(msg{i}.Data, ',');
    if length(val) == 13
        % Convert strings to numbers
        angular_velocity_x = str2double(val{11});
        angular_velocity_y = str2double(val{12});
        calib_gyro_z = split(val{13}, '*');
        angular_velocity_z = str2double(calib_gyro_z{1});
        angular_velocity(i,:) = [angular_velocity_x, angular_velocity_y, angular_velocity_z];
    end
     if  length(val)<=13
         val = zeros(13,1);
     end
end

%plot gyro x,y,z vs time
figure('Name',"Angular Velocity vs Time plot", 'NumberTitle', 'off')
plot(linspace(1,length(msg)/40,length(msg)), angular_velocity(:,1),'-', 'color','r')
hold on
plot(linspace(1,length(msg)/40,length(msg)), angular_velocity(:,2),'-', 'color','g')
hold on
plot(linspace(1,length(msg)/40,length(msg)), angular_velocity(:,3),'-', 'color','b')
grid on
%Labels and titles
xlabel('Time')
ylabel('Angular Velocity')
%legend ('GyroX','GyroY','GyroZ')
title('Angular Velocity vs Time Series')
hold off
saveas(gcf, 'Angular_Velocity_One_Plot.png')



% Allan variance_Angular Velocity_x
Fs = 40;
t0 = 1/Fs;

theta = cumsum(angular_velocity(:,1), 1)*t0;
maxNumM = 100;
L = size(theta, 1);
maxM = 2.^floor(log2(L/2));
m = logspace(log10(1), log10(maxM), maxNumM).';

m = ceil(m); 
m = unique(m); 

tau = m*t0;

[avarFromFunc, tauFromFunc] = allanvar(angular_velocity(:,1), m, Fs);
adevFromFunc = sqrt(avarFromFunc);

% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = -0.5;
logtau = log10(tauFromFunc);
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));

% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);

% Determine the angle random walk coefficient from the line.
logN = slope*log(1) + b;
N_angularVelocity_X = 10^logN

% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = 0.5;
logtau = log10(tauFromFunc);
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the rate random walk coefficient from the line.
logK = slope*log10(3) + b;
K_angularVelocity_X = 10^logK

% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = 0;
logtau = log10(tauFromFunc); %tau = tauFromFunc
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the bias instability coefficient from the line.
scfB = sqrt(2*log(2)/pi);
logB = b - log10(scfB);
B_angularVelocity_X = 10^logB

% Plot
tauN = 1;
lineN = N_angularVelocity_X ./ sqrt(tau);
tauK = 3;
lineK = K_angularVelocity_X .* sqrt(tau/3);
tauB = tau(i);
lineB = B_angularVelocity_X * scfB * ones(size(tau));
% Plot of Allan Deviantion of Gyro in x direction
figure('Name',"Allan Deviation of Gyro X (With Noise Parameters)", 'NumberTitle', 'off')
loglog(tauFromFunc, adevFromFunc, tauFromFunc, lineN, '--', tauN, N_angularVelocity_X, 'o', tauFromFunc, lineK, '--', tauK, K_angularVelocity_X, 'o', tau, lineB, '--', tauB, scfB*B_angularVelocity_X, 'o');
legend('\sigma (rad/s)', '\sigma_N ((rad/s)sqrt{Hz})', '','\sigma_K ((rad/s)sqrt{Hz})','','\sigma_B (rad/s)')
title('Allan Deviation with Noise Parameters of Gyro X')
xlabel('\tau')
ylabel('\sigma(\tau)')
text(tauN, N_angularVelocity_X, 'N')
text(tauK, K_angularVelocity_X, 'K')
text(tauB, scfB*B_angularVelocity_X, '0.664B')
grid on
axis equal
saveas(gcf, 'LocationB_angular_velocity_X.png')





% Allan variance - gyro y
Fs = 40; 
t0 = 1/Fs;
theta = cumsum(angular_velocity(:,2), 1)*t0;
maxNumM = 100;
L = size(theta, 1);
maxM = 2.^floor(log2(L/2));
m = logspace(log10(1), log10(maxM), maxNumM).';
m = ceil(m); 
m = unique(m); 
tau = m*t0;
[avarFromFunc, tauFromFunc] = allanvar(angular_velocity(:,2), m, Fs);
adevFromFunc = sqrt(avarFromFunc);

% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = -0.5;
logtau = log10(tauFromFunc);
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.

b = logadev(i) - slope*logtau(i);
% Determine the angle random walk coefficient from the line.
logN = slope*log(1) + b;
N_angularVelocity_Y = 10^logN

% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = 0.5;
logtau = log10(tauFromFunc);
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the rate random walk coefficient from the line.
logK = slope*log10(3) + b;
K_angularVelocity_Y = 10^logK

% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = 0;
logtau = log10(tauFromFunc); %tau = tauFromFunc
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the bias instability coefficient from the line.
scfB = sqrt(2*log(2)/pi);
logB = b - log10(scfB);
B_angularVelocity_Y = 10^logB

% Plot
tauN = 1;
lineN = N_angularVelocity_Y ./ sqrt(tau);
tauK = 3;
lineK = K_angularVelocity_Y .* sqrt(tau/3);
tauB = tau(i);
lineB = B_angularVelocity_Y * scfB * ones(size(tau));
% Plot of Allan Deviantion of Gyro in y direction
figure('Name',"Allan Deviation of Gyro in Y Direction along with Noise Parameters", 'NumberTitle', 'off')
loglog(tauFromFunc, adevFromFunc, tauFromFunc, lineN, '--', tauN, N_angularVelocity_Y, 'o', tauFromFunc, lineK, '--', tauK, K_angularVelocity_Y, 'o', tau, lineB, '--', tauB, scfB*B_angularVelocity_Y, 'o');
legend('\sigma (rad/s)', '\sigma_N ((rad/s)sqrt{Hz})', '','\sigma_K ((rad/s)sqrt{Hz})','','\sigma_B (rad/s)')
title('Allan Deviation with Noise Parameters of Gyro Y')
xlabel('\tau')
ylabel('\sigma(\tau)')
text(tauN, N_angularVelocity_Y, 'N')
text(tauK, K_angularVelocity_Y, 'K')
text(tauB, scfB*B_angularVelocity_Y, '0.664B')
grid on
axis equal
saveas(gcf, 'LocationB_angular_velocity_Y.png')


% Allan variance - gyro z
Fs = 40;
t0 = 1/Fs;
theta = cumsum(angular_velocity(:,3), 1)*t0;

maxNumM = 100;
L = size(theta, 1);
maxM = 2.^floor(log2(L/2));
m = logspace(log10(1), log10(maxM), maxNumM).';
m = ceil(m); 
m = unique(m); 
tau = m*t0;
angular_velocity(402372,3) = double(0.0);
angular_velocity(402371,3) = double(0.0);

[avarFromFunc, tauFromFunc] = allanvar(angular_velocity(:,3), m, Fs);
adevFromFunc = sqrt(avarFromFunc);
% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = -0.5;
logtau = log10(tauFromFunc);
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the angle random walk coefficient from the line.
logN = slope*log(1) + b;
N_angularVelocity_Z = 10^logN

% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = 0.5;
logtau = log10(tauFromFunc);
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the rate random walk coefficient from the line.
logK = slope*log10(3) + b;
K_angularVelocity_Z = 10^logK

% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = 0;
logtau = log10(tauFromFunc); %tau = tauFromFunc
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the bias instability coefficient from the line.
scfB = sqrt(2*log(2)/pi);
logB = b - log10(scfB);
B_angularVelocity_Z = 10^logB

% Plot
tauN = 1;
lineN = N_angularVelocity_Z ./ sqrt(tau);
tauK = 3;
lineK = K_angularVelocity_Z .* sqrt(tau/3);
tauB = tau(i);
lineB = B_angularVelocity_Z * scfB * ones(size(tau));
% Plot of Allan Deviantion of gyro in Z direction
figure('Name',"Allan Deviation of Gyro in Z Direction along with Noise Parameters", 'NumberTitle', 'off')
loglog(tauFromFunc, adevFromFunc, tauFromFunc, lineN, '--', tauN, N_angularVelocity_Z, 'o', tauFromFunc, lineK, '--', tauK, K_angularVelocity_Z, 'o', tau, lineB, '--', tauB, scfB*B_angularVelocity_Z, 'o');
legend('\sigma (rad/s)', '\sigma_N ((rad/s)sqrt{Hz})', '','\sigma_K ((rad/s)sqrt{Hz})','','\sigma_B (rad/s)')
title('Allan Deviation with Noise Parameters of Gyro Z')
xlabel('\tau')
ylabel('\sigma(\tau)')
text(tauN, N_angularVelocity_Z, 'N')
text(tauK, K_angularVelocity_Z, 'K')
text(tauB, scfB*B_angularVelocity_Z, '0.664B')
grid on
axis equal
saveas(gcf, 'LocationB_angular_velocity_Z.png')
grid on