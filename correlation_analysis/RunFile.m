close all;
clear all;
clc;

% The basic parameters from the experiment and also the coarse graining 
% parameters

% Input basicName = '../xyuv_components/xyuv_'; % Basic prefix for the images
Input.basicName = '../xyuv(64X32)DCC-interpol-2014-03-07/xyuv_';
%Input.basicName = '../PIV_Files_For_Testing_Fine/xyuv_';
Input.format = '.txt'; % format of the image '.txt' '.dat'

% Details of grid size
Input.dx = 35.76; % grid size in 'x' in micrometers
Input.dy = 35.76; % grid size in 'y' in micrometers
Input.dt = 10; % time difference between two frames in minutes


% Input frame details
Input.FrameBegin = 1; % begin from this frame
Input.NFrames = 288; % total number of frames

% Coarse Graining Data

Input.qRBins = 60; % bins for wavenumber q
Input.wTBins = 60; % bins for frequency w
Input.RBins = 60; % bins for r
Input.TBin = 60; % bins for T

% Generate the correlation data with the input
% FlowData has following components
% FlowData is also saved as a file 'FlowData.mat'
% .
%%
FlowData = GenerateCorrelation(Input);

VelCorr = FlowData.VelCorr;
DivCorr = FlowData.DivCorr;
psdVel = FlowData.psdVel;
psdDiv = FlowData.psdDiv;
XMesh = FlowData.XMesh;
YMesh = FlowData.YMesh;
TMesh = FlowData.TMesh;
qXMesh = FlowData.qXMesh;
qYMesh = FlowData.qYMesh;
wTMesh = FlowData.wTMesh;
VelX = FlowData.VelX;
VelY = FlowData.VelY;
Div = FlowData.Div;

%save('./FlowData_Coarse.mat', 'FlowData');

[Ny, Nx, Nt] = size(qXMesh);

clear FlowData;

%%
% We now Analyze the data

% Write Interpolation function for Corr using 2-D spline

DivCtxy = @(x, y, t)interp2(XMesh(:,:,1), YMesh(:,:,1), DivCorr(:,:,t), x, y, 'spline');
VelCtxy = @(x, y, t)interp2(XMesh(:,:,1), YMesh(:,:,1), VelCorr(:,:,t), x, y, 'spline');

DivCtrth = @(r, theta, t) DivCtxy(r*cos(theta), r*sin(theta), t);
VelCtrth = @(r, theta, t) VelCtxy(r*cos(theta), r*sin(theta), t);

DivCtr = @(r, t) 1/(2*pi)*integral(@(theta) DivCtrth(r, theta, t),0, 2*pi, 'AbsTol', 1e-6, 'RelTol', 1e-5);
VelCtr = @(r, t) 1/(2*pi)*integral(@(theta) VelCtrth(r, theta, t),0, 2*pi, 'AbsTol', 1e-6, 'RelTol', 1e-5);

% We use the above interpolation function using 2-D spline
% We now obtain the overall correlation function C(r, t) for Div and Vel


% We first obtain spatial time function

TimeZero = (Nt - 1)/2 + 1; % zero time 

R = sqrt(XMesh(:,:,1).^2 + YMesh(:,:,1).^2);
R = R(:);

NRPoints = 72; % Number of Points to Plot
RBin = linspace(0, max(R), NRPoints); 

RDivCorr = RBin*0;
RVelCorr = RBin*0;

clear R;

for i = 1:length(RBin)
    RDivCorr(i) = DivCtr(RBin(i), TimeZero);
    RVelCorr(i) = VelCtr(RBin(i), TimeZero);
end

RDivCorr = RDivCorr/RDivCorr(1);
RVelCorr = RVelCorr/RVelCorr(1);

NPlot = ceil(NRPoints*0.8);
RDivCorr = RDivCorr(1:NPlot);
RVelCorr = RVelCorr(1:NPlot);
RBin = RBin(1:NPlot);



% Plot Spatial Power Spectrum for Divergence
 figure(1);
 plot(RBin, RDivCorr,'Color','r','Marker','x','LineWidth',2); set(gca,'fontsize',25);
 xlabel('r [\mum]');
 ylabel('Divergence Correlation');
 %set(gca,'XTick',0:0.1:max(wT))
 %saveas(gca, 'GridWaveLengthDiv.eps', 'epsc');
 saveas(gca, 'DivCorr-length.fig')
 saveas(gca, 'DivCorr-length','tif');
 saveas(gca, 'DivCorr-length','svg');
 xdata = RBin; % redefining for minima finding and intercept finding
 ydata = RDivCorr; % redefining for minima finding and intercept finding
 g = @(x)interp1(xdata,ydata,x,'spline');
 fzero(g,250);
 txt=fopen('DivCorr-length (intercept).txt','a+');
 fprintf(txt, ['\r\n','DivCorr-length (intercept) in micrometers = ']);
 fprintf(txt, '%f \t', fzero(g,250)');
 [ydatamax,imax,ydatamin,imin] = extrema(ydata); % To find the minima (using extrema function)
 rDivminima = xdata (imin);
 txt=fopen('DivCorr-length (minima).txt','a+');
 fprintf(txt, ['\r\n','DivCorr-length (minima) in micrometers = ']);
 fprintf(txt, '%f \t', rDivminima');
 

 
 
% Plot Spatial Power Spectrum for Velocity
 figure(2);
 plot(RBin, RVelCorr, 'Color','k','Marker','x','LineWidth',2); set(gca,'fontsize',25);
 xlabel('r [\mum]');
 ylabel('Velocity Correlation');
 %set(gca,'XTick',0:0.1:max(wT))
 %saveas(gca, 'GridWaveLengthDiv.eps', 'epsc');
 saveas(gca, 'VelCorr-length.fig')
 saveas(gca, 'VelCorr-length','tif');
 saveas(gca, 'VelCorr-length','svg');
 xdata = RBin; % redefining for minima finding and intercept finding
 ydata = RVelCorr; % redefining for minima finding and intercept finding
  g = @(x)interp1(xdata,ydata,x,'spline');
 fzero(g,250);
 txt=fopen('VelCorr-length (intercept).txt','a+');
 fprintf(txt, ['\r\n','VelCorr-length (intercept) in micrometers = ']);
 fprintf(txt, '%f \t', fzero(g,250)'); 
 [ydatamax,imax,ydatamin,imin] = extrema(ydata); % To find the minima (using extrema function)
 rVelminima = xdata (imin);
 txt=fopen('VelCorr-length.txt','a+');
 fprintf(txt, ['\r\n','VelCorr-length (minima) in micrometers = ']);
 fprintf(txt, '%f \t', rVelminima');
 
 
 
 
%%
% We now obtain temporal correlation
 
TimeZero = (Nt - 1)/2 + 1; % zero time 
TBin = [0:TimeZero-1]*Input.dt/60; % time in hours
% Time correlation
TDivCorr = TBin;
TVelCorr = TBin;


TDivCorr = TBin*0;
TVelCorr = TBin*0;

for i = 1:length(TBin)
    TDivCorr(i) = DivCtr(0, TimeZero+i-1);
    TVelCorr(i) = VelCtr(0, TimeZero+i-1);
end
TDivCorr = TDivCorr/TDivCorr(1);
TVelCorr = TVelCorr/TVelCorr(1);

NPoints = length(TBin);

NPlot = ceil(NRPoints*0.8);
TDivCorr = TDivCorr(1:NPlot);
TVelCorr = TVelCorr(1:NPlot);
TBin = TBin(1:NPlot);


 figure(3);
 plot(TBin, TDivCorr, 'Color','r','Marker','o','LineWidth',2); set(gca,'fontsize',25);
 xlabel('Time [h]');
 ylabel('Divergence Correlation');
 %set(gca,'XTick',0:0.1:max(wT))
 %saveas(gca, 'GridWaveLengthDiv.eps', 'epsc');
 saveas(gca, 'DivCorr-period.fig')
 saveas(gca, 'DivCorr-period','tif');
 saveas(gca, 'DivCorr-period','svg');
 xdata = TBin; % redefining for minima finding and intercept finding
 ydata = TDivCorr; % redefining for minima finding and intercept finding
 [ydatamax,imax,ydatamin,imin] = extrema(ydata); % To find the minima (using extrema function)
 %tDivminima = xdata (imin);
 Divergenceperiod =  xdata (imin)*2; % Twice of the obtained minima
 txt=fopen('DivCorr-period.txt','a+');
 %fprintf(txt, ['\r\n','DivCorr-minima in hours = ']);
 %fprintf(txt, '%f \t', tDivminima');
 fprintf(txt, ['\r\n','Divergence period in hours = ']); % Twice of the obtained minima
 fprintf(txt, '%f \t', Divergenceperiod');
 
 figure(4);
 plot(TBin, TVelCorr, 'Color','k','Marker','o','LineWidth',2); set(gca,'fontsize',25);
 xlabel('Time [h]');
 ylabel('Velocity Correlation');
 %set(gca,'XTick',0:0.1:max(wT))
 %saveas(gca, 'GridWaveLengthDiv.eps', 'epsc');
 saveas(gca, 'VelCorr-period.fig')
 saveas(gca, 'VelCorr-period','tif');
 saveas(gca, 'VelCorr-period','svg');
 xdata = TBin; % redefining for minima finding and intercept finding
 ydata = TVelCorr; % redefining for minima finding and intercept finding
 [ydatamax,imax,ydatamin,imin] = extrema(ydata); % To find the minima (using extrema function)
 %tVelminima = xdata (imin);
 Velocityperiod = xdata (imin)*2; % Twice of the obtained minima
 txt=fopen('VelCorr-period.txt','a+');
 %fprintf(txt, ['\r\n','Velcorrminima in hours = ']);
 %fprintf(txt, '%f \t', tVelminima');
 fprintf(txt, ['\r\n','Velocity period in hours = ']); % Twice of the obtained minima
 fprintf(txt, '%f \t', Velocityperiod');
 


%% We now obtain the overall powerspectrum in space and time for vel and div

wT = wTMesh(1,1,:);
wT = wT(:)';
qX = qXMesh(1,:,1);
qY = qYMesh(:, 1, 1)';
qR = sqrt(qXMesh(:,:,1).^2 + qYMesh(:,:,1).^2);
qR = qR(:);
qR = sort(union(qR, qR));


% Write Interpolation function for PSD
pDivwqxqy = @(qx, qy, w)interp2(qXMesh(:,:,1), qYMesh(:,:,1), psdDiv(:,:,w), qx, qy, 'spline');
pVelwqxqy = @(qx, qy, w)interp2(qXMesh(:,:,1), qYMesh(:,:,1), psdVel(:,:,w), qx, qy, 'spline');

pDivwqrqt = @(qr, qt, w) pDivwqxqy(qr*cos(qt), qr*sin(qt), w);
pVelwqrqt = @(qr, qt, w) pVelwqxqy(qr*cos(qt),qr*sin(qt), w);

pDivwqr = @(qr, w) 1/(2*pi)*integral(@(qt) pDivwqrqt(qr, qt, w),0, 2*pi, 'AbsTol', 1e-6, 'RelTol', 1e-5);
pVelwqr = @(qr, w) 1/(2*pi)*integral(@(qt) pVelwqrqt(qr, qt, w),0, 2*pi, 'AbsTol', 1e-6, 'RelTol', 1e-5);


%% Get Total Power spectrum over all wavelengths w.r.t time

 FreqPDiv = zeros(length(wT), 1);
 SpaceDiv = zeros(length(qY), length(qX));
 
 FreqPVel = zeros(length(wT), 1);
 SpaceVel = zeros(length(qY), length(qX));
 
 
 
 for i = 1:length(wT)
     divFreq = psdDiv(:,:,i);
     velFreq = psdVel(:, :, i);
     FreqPDiv(i) = sum(divFreq(1:end));
     FreqPVel(i) = sum(velFreq(1:end));
     %FreqPDiv(i) = pDivwqr(qR(500), i);
     %FreqPVel(i) = pVelwqr(qR(500), i);
     SpaceDiv = SpaceDiv + psdDiv(:,:,i);
     SpaceVel = SpaceVel + psdVel(:, :, i);
     
 end
 SpaceDiv = SpaceDiv/length(wT);
 SpaceVel = SpaceVel/length(wT);
 
 NTotal = length(wT);
 nwZero = (NTotal-1)/2 + 1;
 
 figure(5);
 hold on
 semilogy(wT(nwZero:end), FreqPDiv(nwZero:end), 'color', 'r', 'Marker','o','LineWidth',2); set(gca,'fontsize',25);
 xlabel('Frequency [1/h]');
 ylabel('Mean Divergence PSD'); % Mean over all wavelengths
 %set(gca,'XTick',0:0.1:max(wT))
 %saveas(gca, 'GridWaveLengthDiv.eps', 'epsc');
 saveas(gca, 'DivPSD-period.fig')
 saveas(gca, 'DivPSD-period','tif');
 saveas(gca, 'DivPSD-period','svg');
 [val t] = max(FreqPDiv);
 DivPSDperiod = 1/abs(wT(t)) % time period of the possible oscillation in hours
 txt=fopen('DivPSD-period.txt','w');
 fprintf(txt, ['Period in hours = ' '\n']);
 fprintf(txt, '%f \n', DivPSDperiod');
  
 
 
 figure(6);
 hold on
 semilogy(wT(nwZero:end), FreqPVel(nwZero:end), 'color', 'k', 'Marker','o','LineWidth',2); set(gca,'fontsize',25);
 xlabel('Frequency [1/h]');
 ylabel('Mean Velocity PSD'); % Mean over all wavelengths
 %set(gca,'XTick',0:0.1:max(wT))
 %saveas(gca, 'GridWaveLengthVel.eps', 'epsc');
 saveas(gca, 'VelPSD-period.fig')
 saveas(gca, 'VelPSD-period','tif');
 saveas(gca, 'VelPSD-period','svg');
 [val t] = max(FreqPVel);
 VelPSDperiod = 1/abs(wT(t)) % time period of the possible oscillation in hours
 txt=fopen('VelPSD-period.txt','w');
 fprintf(txt, ['Period in hours = ' '\n']);
 fprintf(txt, '%f \n', VelPSDperiod');
 
%%


pDivwqxqy = @(qx, qy)interp2(qXMesh(:,:,1), qYMesh(:,:,1), SpaceDiv(:,:), qx, qy, 'spline');
pVelwqxqy = @(qx, qy)interp2(qXMesh(:,:,1), qYMesh(:,:,1), SpaceVel(:,:), qx, qy, 'spline');

pDivwqrqt = @(qr, qt) pDivwqxqy(qr*cos(qt), qr*sin(qt));
pVelwqrqt = @(qr, qt) pVelwqxqy(qr*cos(qt),qr*sin(qt));

pDivwqr = @(qr) 1/(2*pi)*integral(@(qt) pDivwqrqt(qr, qt),0, 2*pi, 'AbsTol', 1e-6, 'RelTol', 1e-5);
pVelwqr = @(qr) 1/(2*pi)*integral(@(qt) pVelwqrqt(qr, qt),0, 2*pi, 'AbsTol', 1e-6, 'RelTol', 1e-5);

NPoints = 150;
qRBin = linspace(0, max(qR), NPoints);

psdDivQ = zeros(length(qRBin),1);
psdVelQ = zeros(length(qRBin),1);

for i = 1:length(qRBin)
    psdDivQ(i) = pDivwqr(qRBin(i));
    psdVelQ(i) = pVelwqr(qRBin(i));
end



figure(7)
 semilogy(qRBin(1:100), psdVelQ(1:100), 'color', 'k', 'Marker','x','LineWidth',2); set(gca,'fontsize',25);
 xlabel('Wavelength [1/\mum]');
 ylabel('Mean Velocity PSD'); % Mean over all frequencies
 %set(gca,'XTick',0:0.1:max(wT))
 %saveas(gca, 'GridWaveLengthVel.eps', 'epsc');
 saveas(gca, 'VelPSD-wavelength.fig')
 saveas(gca, 'VelPSD-wavelength','tif');
 saveas(gca, 'VelPSD-wavelength','svg');
 [val w] = max(psdVelQ);
 VelPSDwavelength = 1/abs(qRBin(w)) % Wavelength of the pulsation in micrometers 
 txt=fopen('VelPSD-wavelength.txt','w');
 fprintf(txt, ['Wavelength (Vel) in micrometers = ' '\n']);
 fprintf(txt, '%f \n', VelPSDwavelength');

 
 figure(8)
 semilogy(qRBin(1:100), psdDivQ(1:100), 'color', 'r', 'Marker','x','LineWidth',2); set(gca,'fontsize',25);
 xlabel('Wavelength [1/\mum]');
 ylabel('Mean Divergence PSD'); % Mean over all frequencies
 %set(gca,'XTick',0:0.1:max(wT))
 %saveas(gca, 'GridWaveLengthVel.eps', 'epsc');
 saveas(gca, 'DivPSD-wavelength.fig')
 saveas(gca, 'DivPSD-wavelength','tif');
 saveas(gca, 'DivPSD-wavelength','svg');
 [val w] = max(psdDivQ);
 DivPSDwavelength = 1/abs(qRBin(w)) % Wavelength of the pulsation in micrometers 
 txt=fopen('DivPSD-wavelength.txt','w');
 fprintf(txt, ['Wavelength (Div) in micrometers = ' '\n']);
 fprintf(txt, '%f \n', DivPSDwavelength');

 close all;