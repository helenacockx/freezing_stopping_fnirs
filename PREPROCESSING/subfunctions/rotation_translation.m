function [rot_tra] = rotation_translation(data, run, ID, varargin)
% ROTATION_TRANSLATION
% This function determines the translation vector and rotation angle
% to reframe the motion data into the new coordinate system based on the
% trajectory of the center of mass
%
% Use as
%   [data_reframed, rot_angle, tra_vec] = rotation_translation(data, run, ID, varargin);
%
% INPUT:
%       data         = motion data in fieldtrip data structure
%       run          = run number of which the data was given
%       ID           = ID of which the data was given 
% Additional options can be specified in key-value pairs and can be:
%       'visualize'    true or false (default = true)
%
% OUTPUT
%       rot_tra       = structure containing info about the translation and rotation, containing the following fields:
%         tra_vec     = translation vector [x y z] that was used to move the center
%                       of the trajectory in the x, y and z direction to the midpoint.
%         rot_angle   = angle that was used to rotate the motion data  
%         corr_door   = correction factor that was used to translate the x
%                       position such that zerocrossing was around the door
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  get the options
vis=ft_getopt(varargin, 'visualize', true);

%% calculate translation vector based on the Center of Mass (COM) trajectory
% define the trajectory of the COM
COM.pos.X = data.trial{1}(find(startsWith(data.label, 'seg_COM_centerOfMass_X')), :);
COM.pos.Y = data.trial{1}(find(startsWith(data.label, 'seg_COM_centerOfMass_Y')), :);
COM.pos.Z = data.trial{1}(find(startsWith(data.label, 'seg_COM_centerOfMass_Z')), :);

if vis
    figure; axis vis3d; axis on; grid on; hold on; view(0,90);
    plot3(COM.pos.X, COM.pos.Y, COM.pos.Z, 'b.');
    title(sprintf('trajectory of center of mass in absolute positions (sub-%s run-%02d)', ID, run))
end

% find the translation vector to move the middle of the COM trajectory to
% the (0,0) coordinate
tra_vec=[-(max(COM.pos.X)+min(COM.pos.X))/2 -(max(COM.pos.Y)+min(COM.pos.Y))/2 0];

% translate COM trajectory to check
COM_tra.pos.X=COM.pos.X + tra_vec(1);
COM_tra.pos.Y=COM.pos.Y + tra_vec(2);
COM_tra.pos.Z=COM.pos.Z + tra_vec(3);

if vis
    hold on; plot3(COM_tra.pos.X, COM_tra.pos.Y, COM_tra.pos.Z, 'm.');
end

%% calculate the rotation angle based on the center of mass (COM) trajectory
% select the point with the highest x coordinate and define rotation angle
if (contains(ID, {'HC35', 'HC76', 'PD25'}) & run==3) | (contains(ID, 'HC34') & run==2) % exception
    [m,i]=max(COM_tra.pos.X);
    x_max=COM_tra.pos.X(1,i);
    y_max=COM_tra.pos.Y(1,i);
else
  [pks, locs]=findpeaks(COM_tra.pos.X, 'MinPeakHeight', 2, 'MinPeakDistance', 30*data.fsample);
  x_max=median(pks);
  y_max=median(COM_tra.pos.Y(locs));
end
a1 = atand(y_max/x_max);

% same for the lowest x coordinate
if contains(ID, {'HC35', 'HC76', 'PD25'}) & run==3 % exception
  [m,i]=min(COM_tra.pos.X);
  x_min=COM_tra.pos.X(1,i);
  y_min=COM_tra.pos.Y(1,i);
elseif contains(ID, 'HC34') & run == 2
    [m,i]=min(COM_tra.pos.Y);
  x_min=COM_tra.pos.X(1,i);
  y_min=COM_tra.pos.Y(1,i);
else
  [pks, locs]=findpeaks(-COM_tra.pos.X, 'MinPeakHeight', 2, 'MinPeakDistance', 30*data.fsample);
  x_min=median(-pks);
  y_min=median(COM_tra.pos.Y(locs));
end
a2 = atand(y_min/x_min);

% take the mean of both and reverse (we want to rotate in the opposite
% direction
rot_angle=-mean([a1 a2]);
if contains(ID, 'PD10') & run==2 % exception: weird angle at x max
  rot_angle=-a2;
elseif contains(ID, 'HC34') & run==2
  rot_angle=mean([-a1 a2]);
end

% rotate COM trajectory  to check
COM_rot.pos.X = COM_tra.pos.X*cosd(rot_angle)-COM_tra.pos.Y*sind(rot_angle);
COM_rot.pos.Y = COM_tra.pos.X*sind(rot_angle)+COM_tra.pos.Y*cosd(rot_angle);
COM_rot.pos.Z = COM_tra.pos.Z;

if vis
    hold on; plot3(COM_rot.pos.X, COM_rot.pos.Y, COM_rot.pos.Z, 'r.');
    plot(x_max, y_max, 'ko'); plot(x_min, y_min, 'ko')
    legend('original', 'translated', 'translated + rotated')
end

%% calculate the translation in the X axis based on the doorway crossing
stand=find(([0 diff(COM_rot.pos.X)]<0.001 & [0 diff(COM_rot.pos.X)]>-0.001) & (COM_rot.pos.X<5 & COM_rot.pos.X>-5));
stand=stand(find(stand>70*data.fsample)); % do not include the standing period at the start
corr_door=mean(COM_rot.pos.X(stand)); 
if strcmp(ID, 'PD90') & run==4
  corr_door=0.4;
elseif strcmp(ID, 'PD88') & run==4
  corr_door=0;
end
if vis
  figure; 
  plot(COM_rot.pos.X); grid on; title(sprintf('center of mass trajectory along x-axis before doorway correction (sub-%s run-%02d)', ID, run))
  hold on; plot(stand, zeros(1, length(stand)), 'ro')
  plot(corr_door*ones(1,length(COM_rot.pos.X)))
  legend({'COM', 'doorway standing', 'calculated doorway position'})
end
tra_vec(1)=tra_vec(1);

%% save the translation & rotation info in a structure
rot_tra.tra_vec=tra_vec;
rot_tra.rot_angle=rot_angle;
rot_tra.corr_door=corr_door;
    
