% MMI: @MPI-PKS, Dresden, 18th June 2015

function VelData = ExtractVel(FileName, dx, dy)
CompFile = importdata(FileName); % load the complete file with headers

FileType = whos('CompFile');
if (FileType.class == 'struct')
    D = CompFile.data; % Use only the data from the file
else
    D = CompFile;
end

clear CompFile; % clear the huge imported file

%D = load('/Users/mandarinamdar/Dropbox/Research/Guillaume/Raghavan-Winding-Number/PosVel.dat');

xx = D(:,1);% x position
maxX = max(xx); % max value of x
yy = D(:,2); % y position
maxY = max(yy); % top co-ordinate of yy
yy = maxY - yy; % this is to modify the up down symmtery
uu = D(:,3); % x component of velocity
vv = -D(:,4); % y component of velocity
%ss = sqrt(uu.^2 + vv.^2); % speed

% subtracting the average velocities in the x and y direction for a frame
%uu = uu - mean(uu(1:end));
%vv = vv - mean(vv(1:end));
% define a grid as per the experimental data
xG = min(xx):dx:max(xx);
yG = min(yy):dy:max(yy);

[Xgrid, Ygrid] = meshgrid(xG, yG);

clear xG yG; % clear the temporary variables

% reformat all the points on the grid points systematically
% obtain x and y components of velocity and divergence
velX = griddata(xx, yy, uu, Xgrid, Ygrid);
velY = griddata(xx, yy, vv, Xgrid, Ygrid);
div = divergence(Xgrid, Ygrid, velX, velY);

clear xx yy uu vv; % clear more temporary variables

VelData = cell(5,1); % cell object for storing the details
VelData{1} = Xgrid;
VelData{2} = Ygrid;
VelData{3} = velX;
VelData{4} = velY;
VelData{5} = div;

clear Xgrid Ygrid velX velY;

% Saving the quiver plot to an Image
% %q = quiver(VelData{1}, VelData{2}, VelData{3}, VelData{4});
% set(q, 'color', 'black', 'linewidth', 1);
% title('Velocity vectors');
% xlabel('X');
% ylabel('Y');
% set(gca,'FontSize',20)
% set(findall(gcf,'type','text'),'FontSize',20)
% box on;
% xlim([0-dx, maxX + dx]);
% ylim([0-dy, maxY + dy]);
% EpsName = strcat(FileName(1:end-3), 'tif');
% saveas(gca, EpsName,'tif');

end