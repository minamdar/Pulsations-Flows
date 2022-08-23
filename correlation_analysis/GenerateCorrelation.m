
% History of modification
% MMI: IIT Bombay, 19 September 2016
% MMI: IIT Bombay, 10 June 2016

% MMI: IIT Bombay
% Modifications made on 20 May 2016

% MMI, @MPI-PKS, Dresden, 19th June 2015
% variable 't' is integer variable for time
% 'i' is integer variable for y
% 'j' is integer variable for x

% What does this code do:
% Loads images in the directory 
% Finds correlatin function and powerspectrum for divergence and velocity

%%
function [FlowData] = GenerateCorrelation(Input)
% Read the input

close all;

%basicName = '../xyuv_components/xyuv_'; % Basic prefix for the images
basicName = Input.basicName;
format = Input.format; % format of the image '.txt' '.dat'

% Details of grid size
dx = Input.dx; % grid size in 'x'
dy = Input.dy; % grid size in 'y'
dt = Input.dt; % time difference between two frames in minutes


% Creating the X and Y Grid and the corresponding details
FileName = strcat(basicName, sprintf('%04d',1), '.txt')
V = ExtractVel(FileName, dx, dy);
Xgrid = V{1}; % x grid: same for all times
Ygrid = V{2}; % y grid: same for all times
[NGridY, NGridX] = size(Xgrid); % size of the grid in X and Y direction
Lx = dx*(NGridX-1); % size of the domain in X
Ly = dy*(NGridY-1); % size of domain in Y

clear V; % clear the variable

% Loading frame details
FrameBegin = Input.FrameBegin; % begin from this frame
NFrames = Input.NFrames; % total number of frames
FrameEnd = FrameBegin + NFrames; % final frame


% Creating arrays for velocity and divergence
velx = zeros(NGridY, NGridX, NFrames); % Define a 3D array
vely = zeros(NGridY, NGridX, NFrames); % Define a 3D array
div = zeros(NGridY, NGridX, NFrames); % Define a 3D array


for t = 1:NFrames
    FrameNumber = t+FrameBegin-1;
    FileName = strcat(basicName, sprintf('%04d',FrameNumber), '.txt'); % construct the filename
    V = ExtractVel(FileName, dx, dy);
    velx(:,:,t) = V{3}; % x vel at i time
    vely(:,:,t) = V{4}; % y vel at i time
    div(:,:,t) = V{5}; % div at i time  
    % cell variables: i do not know how to manipulate them nicely
%     VelDiv{i, 1} = V{3}; % x-vel for i time
%     VelDiv{i, 2} = V{4}; % y-vel for i time
%     VelDiv{i, 3} = V{5}; % div for i time
end


% *** All the manipulations below are only to ensure that
% *** there has been one-to-one correspondence between the 
% *** the power spectrum and FT of the correlation

%*** We pad the VelX, VelY and Div with zeros in the time dimension
VelX = padarray(velx, [0,0,NFrames-1]);
VelY = padarray(vely, [0,0,NFrames-1]);
Div = padarray(div, [0,0,NFrames-1]);


%*** We obtain the size of the basic array
sizeV = size(VelX);


% We obtain the correlation functions
% We will first get C(x,y,t) for divergence and velocity
% Note that the if Nx, Ny, Nz is the size of the original array
% the size of the correlation function array is 2*Nx-1, 2*Ny-1, 2*Nz-1

VelCorr = zeros(2*NGridY-1, 2*NGridX-1, 2*NFrames-1);
DivCorr = zeros(2*NGridY-1, 2*NGridX-1, 2*NFrames-1);

%***Here we obtain the correlation function 
% ** Since Matlab only allows maximum of xcorr2
% ** we need to write one loop over time 't'

for t = -(NFrames-1):1:(NFrames-1)
    vcorrtemp = zeros(2*NGridY-1, 2*NGridX-1);
    divcorrtemp = zeros(2*NGridY-1, 2*NGridX-1);
    for j = 1:NFrames
        VXTemp = VelX(:,:,t+NFrames+j-1);
        vxtemp = velx(:,:,j);
        
        VYTemp = VelY(:,:,t+NFrames+j-1);
        vytemp = vely(:,:,j);
        
        DivTemp = Div(:,:,t+NFrames+j-1);
        divtemp = div(:,:,j);
        
       % aTemp = a(:,:,j);
        vcorrtemp = vcorrtemp + xcorr2(VXTemp, vxtemp) + ...
                    xcorr2(VYTemp, vytemp);
        
        divcorrtemp = divcorrtemp + xcorr2(DivTemp, divtemp);
        
    end
    %size(temp)
    %temp = circshift(temp, [-(Ny-1), -(Nx-1)]);
    disp(['  ' num2str(t) ' of ' num2str(NFrames) ' lines...']);  
    DivCorr(:, :, t+NFrames) = divcorrtemp;
    VelCorr(:, :, t+NFrames) = vcorrtemp;
end


%**** Here we do check between the PSD and FT of Correlation
% ** We define a array DivCorrA which where we shift towards top left
% ** by an amount of N-1
% ** We also pad 'array' with zeros to finally match it with the dimension 
% ** of Correlation 


%** Here is for divergence correlation
DivCorrA = DivCorr;
DivCorrA = circshift(DivCorrA, [-(NGridY-1), -(NGridX-1), -(NFrames-1)]);
DivA = padarray(div, [ceil((NGridY-1)/2), ceil((NGridX-1)/2), ceil((NFrames-1)/2)]);


%** Here is for velocity correlation
VelCorrA = VelCorr;
VelCorrA = circshift(VelCorrA, [-(NGridY-1), -(NGridX-1), -(NFrames-1)]);
VelAX = padarray(velx, [ceil((NGridY-1)/2), ceil((NGridX-1)/2), ceil((NFrames-1)/2)]);
VelAY = padarray(vely, [ceil((NGridY-1)/2), ceil((NGridX-1)/2), ceil((NFrames-1)/2)]);


% ** If any of the Nx, Ny, Nt is even then we have to remove one row of
% ** padding

if(mod(NGridY,2)==0)
    DivA = DivA(1:end-1,:,:);
    VelAX = VelAX(1:end-1,:,:);
    VelAY = VelAY(1:end-1,:,:);
end

if(mod(NGridX,2)==0)
    DivA = DivA(:,1:end-1,:);
    VelAX = VelAX(:,1:end-1,:);
    VelAY = VelAY(:,1:end-1,:);
end

if(mod(NFrames,2)==0)
    DivA = DivA(:,:,1:end-1);
    VelAX = VelAX(:,:,1:end-1);
    VelAY = VelAY(:,:,1:end-1);
end

% *** Here we obtain the power spectrum of the matrix with padded zeros
% ** power spectrum for Divergence is 'psdDiv'
% ** power spectrum for Velocity is 'psdVel'

psdDiv = abs(fftn(DivA)).^2;
psdVel = abs(fftn(VelAX)).^2 + abs(fftn(VelAY)).^2;


% ** We obtain the Fourier Transform for Padded Divergence and Velocity
% ** correlation

fftCorrDiv = abs(fftn(DivCorrA));
fftCorrVel = fftn(VelCorrA);


%%** Checking the correspondence between FT and PSD
diffDiv = psdDiv - fftCorrDiv;

max(abs(diffDiv(:)));
max(psdDiv(:));
max(fftCorrDiv(:));

diffVel = psdVel - fftCorrVel;


% Saving the Grid data both for space-time and q-w space

% Obtain dimensions of the matrices to store
[Ny Nx Nt] = size(DivCorr);

% Generate the wave-number data
qX = [-(Nx-1)/2:(Nx-1)/2]/(Nx*dx);
qY = [-(Ny-1)/2:(Ny-1)/2]/(Ny*dy);
wT = [-(Nt-1)/2:(Nt-1)/2]*60/(Nt*dt);

% Generate the space mesh data
X = [-(Nx-1)/2:(Nx-1)/2]*dx;
Y = [-(Ny-1)/2:(Ny-1)/2]*dy;
T = [-(Nt-1)/2:(Nt-1)/2]*dt/60; % in hours

% Now get meshgrid in the q-w space 
[qXMesh, qYMesh, wTMesh] = meshgrid(qX, qY, wT);

[XMesh, YMesh, TMesh] = meshgrid(X, Y, T);

% Obtain the power-spectrum with zero at the middle and symmetric
% on other sides

psdDiv = circshift(psdDiv, [(Ny-1)/2, (Nx-1)/2, (Nt-1)/2]);
psdVel = circshift(psdVel, [(Ny-1)/2, (Nx-1)/2, (Nt-1)/2]);

% Saving the appropriate files in cell variables. 

FlowData.VelCorr = VelCorr;
FlowData.DivCorr = DivCorr;
FlowData.psdVel = psdVel;
FlowData.psdDiv = psdDiv;
FlowData.XMesh = XMesh;
FlowData.YMesh = YMesh;
FlowData.TMesh = TMesh;
FlowData.qXMesh = qXMesh;
FlowData.qYMesh = qYMesh;
FlowData.wTMesh = wTMesh;
FlowData.VelX = velx;
FlowData.VelY = vely;
FlowData.Div = div;

% saving the flow data
end