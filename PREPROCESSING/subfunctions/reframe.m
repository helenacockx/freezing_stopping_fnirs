function [data_reframed]=reframe(data, chan, reframe, rot_tra)
% REFRAME
% This function uses the translation vector, rotation angle and door
% correction of ROTATION_TRANSLATION to reframe the input data into the new coordinate system
%
% Use as
%   [data_reframed] = reframe(data, chan, rot_tra);
%
% INPUT:
%       data         = motion data in fieldtrip data structure that you
%       would like to reframe to the new coordinate system
%       chan         = channel names of the data you would like to reframe
%       (without appendix X/Y/Z)
%       reframe      = 'rot+tra' or 'rot' for rotation + translation
%       (position data) or rotation only, respectively (velocity,acceleration,...)
%       rot_tra      = rotation and translation parameters outputted from
%       ROTATION_TRANSLATION containing the following fields:
%         tra_vector   = translation vector [x y z] that was outputted from
%         ROTATION_TRANSLATION
%         rot_angle    = angle that was outputted from ROTATION_TRANSLATION
%         corr_door    = correction factor outputted from ROTATION_TRANSLATION
%         that was used to translate the x position such that zerocrossing was
%         around the door
%       run          = run number of which the data was given
%       ID           = ID of which the data was given 
%
% OUTPUT:
%       data_reframed = position data in fieldtrip format adapted to the
%       new coordinate system
% 
% DEPENDENCIES:/

%% get the options
tra_vec=rot_tra.tra_vec;
rot_angle=rot_tra.rot_angle;
corr_door=rot_tra.corr_door;

%% select the data of the given channel
cfg=[];
cfg.channel=find(contains(data.label, chan));
data_sel=ft_selectdata(cfg,data);

%% rotate and translate the data
data_reframed=data_sel;
for i=1:3:length(data_sel.label)-2
  if strcmp(reframe, 'rot+tra')
    fprintf('rotating and translating %s... \n', data_sel.label{i:i+2})
    data_tra_X=data_sel.trial{1}(i,:)+tra_vec(1);
    data_tra_Y=data_sel.trial{1}(i+1,:)+tra_vec(2);
    data_tra_Z=data_sel.trial{1}(i+2,:)+tra_vec(3);
    data_reframed.trial{1}(i,:)=data_tra_X*cosd(rot_angle)-data_tra_Y*sind(rot_angle)-corr_door; % -corr_door
    data_reframed.trial{1}(i+1,:)=data_tra_X*sind(rot_angle)+data_tra_Y*cosd(rot_angle);
    data_reframed.trial{1}(i+2,:)=data_tra_Z;
  elseif strcmp(reframe, 'rot')
    data_tra_X=data_sel.trial{1}(i,:);
    data_tra_Y=data_sel.trial{1}(i+1,:);
    data_tra_Z=data_sel.trial{1}(i+2,:);
    data_reframed.trial{1}(i,:)=data_tra_X*cosd(rot_angle)-data_tra_Y*sind(rot_angle);
    data_reframed.trial{1}(i+1,:)=data_tra_X*sind(rot_angle)+data_tra_Y*cosd(rot_angle);
    data_reframed.trial{1}(i+2,:)=data_tra_Z;
  end
end
