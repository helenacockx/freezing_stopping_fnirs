function [run]=load_run_info(run, ID)
%% LOAD_RUN_INFO
% Loads information about the runs for each participant. E.g. if the runs
% are complete or not & what the valid window for this run is.
%
% INPUT:
% - run: matlab structure containing the raw nirs & motion data of each run
% - ID: participant-ID of the run that is currently loaded
%
% OUTPUT:
% - run: same as input, but info added about each run in the "info" field
% 
% DEPENDENCIES: none
%
%% 
% add run info
    switch ID
      case {'PD06', 'PD10', 'PD12','PD22', 'PD35', 'PD46', 'PD48', 'PD50','PD57', 'PD62', 'PD63', 'PD76', 'PD77', 'PD17', 'PD45',...
          'HC01', 'HC04', 'HC06', 'HC13', 'HC19', 'HC28', 'HC33', 'HC35', 'HC41', 'HC42', 'HC60', 'HC66', 'HC67', 'HC68', 'HC76', 'HC90', 'HC91'}
        for i=1:length(run)
          if isempty(run(i).data_motion) | isempty(run(i).data_nirs)
            run(i).info.complete = 'no';
            run(i).info.validwindow = [0 0];
            run(i).info.remark = 'nirs/motion data is missing';
          else
            run(i).info.complete='yes';
            run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
          end
        end
      case 'PD11'
        run(1).info.complete='no';
        run(1).info.validwindow=[0 104];
        run(1).info.remark='sensor left upper leg loose in last seconds'; 
        for i=[2:5]
          run(i).info.complete='yes';
          run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
        end
      case 'PD15'
        run(4).info.complete='no';
        run(4).info.validwindow=[0 418];
        run(4).info.remark='xsens signal got lost';
        for i=1:3
          run(i).info.complete='yes';
          run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
        end
      case 'PD16'
        run(1).info.complete = 'yes';
        run(1).info.validwindow=[0 run(1).data_nirs.data_raw.time{1}(end)];
        run(2).info.complete = 'yes';
        run(2).info.validwindow = [0 run(2).data_nirs.data_raw.time{1}(end)];
      case 'PD25'
        run(1).info.complete = 'no';
        run(1).info.validwindow = [0 180];
        run(1).info.remark='xsens signal got lost';
        run(3).info.complete = 'no';
        run(3).info.validwindow = [0 260];
        run(3).info.remark='xsens signal got lost';
        for i=[2 4:6]
          run(i).info.complete = 'yes';
          run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
        end
      case 'PD31'
        run(2).info.complete='no';
        run(2).info.validwindow=[0 148];
        run(2).info.remark='sensor right upper leg loose in last seconds';
        for i=[1 3:5]
          run(i).info.complete='yes';
          run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
        end
      case 'PD57'
        run(4).info.complete = 'no';
        run(4).info.validwindow  = [0 504];
        run(4).info.remark = 'nirs data went missing at the end'; 
        for i=1:3
          run(i).info.complete = 'yes';
          run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
        end
      case 'PD61'
        run(1).info.validwindow=[0 208];
        run(1).info.complete='no';
        run(1).info.remark='sensor left upper leg loose during whole run';
        for i=2:5
          run(i).info.complete='yes';
          run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
        end
      case 'PD88'
        for i=1:3
          run(i).info.complete='yes';
          run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
        end
        run(4).info.validwindow=[0 150];
        run(4).info.complete='no';
        run(4).info.remark='ended run because participant needed to go to the toilet';
      case 'PD89'
        run(1).info.validwindow=[0 352];
        run(1).info.complete= 'no';
        run(1).info.remark = 'sensor left foot fell off';
        for i=2:5
          run(i).info.complete='yes';
          run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
        end
      case 'PD90'
        run(1).info.validwindow=[0 0];
        run(1).info.complete='no';
        run(1).info.remark= 'bad callibration';
        
        run(4).info.validwindow=[0 125];
        run(4).info.complete='no';
        run(4).info.remark='almost fell at the end';
        cfg=[];
        cfg.latency=run(4).info.validwindow;
        
        run(5).info.validwindow=[59 run(5).data_motion.data_raw.time{1}(end)];
        run(5).info.complete='no';
        run(5).info.remark='was not standing still at the beginning of the run';
        for i=[2 3 6]
          run(i).info.complete='yes';
          run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
        end
      case 'PD96'
        run(2).info.validwindow=[0 0];
        run(2).info.complete='no';
        run(2).info.remark='xsens signal got lost';
        run(2).data_motion=[];
        for i=[1 3:5]
          run(i).info.complete='yes';
          run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
        end
      case 'HC02'
        run(1).info.validwindow = [0 400];
        run(1).info.complete = 'no';
        run(1).info.remark = 'xsens signal went drifting at the end';
        for i=2:4
          run(i).info.complete='yes';
          run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
        end
      case 'HC21'
        run(1).info.validwindow = [0 210];
        run(1).info.complete = 'no';
        run(1).info.remark = 'xsens signal went drifting';
        for i=2:5
          run(i).info.complete='yes';
          run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
        end          
      case 'HC29'
        run(4).info.validwindow = [0 360];
        run(4).info.complete = 'no';
         run(4).info.remark = 'xsens signal went drifting';
        for i=1:3
          run(i).info.complete='yes';
          run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
        end            
      case 'HC34'
        run(4).info.validwindow = [0 170];
        run(4).info.complete = 'no';
        run(4).info.remark = 'xsens signal went drifting';
        for i=1:3
          run(i).info.complete='yes';
          run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
        end
      case 'HC57'
        run(1).info.validwindow = [0 162];
        run(1).info.complete = 'no';
        run(1).info.remark = 'xsens sensor fell off';
        for i=2:6
          if isempty(run(i).data_motion) | isempty(run(i).data_nirs)
            run(i).info.complete = 'no';
            run(i).info.validwindow = [0 0];
            run(i).info.remark = 'nirs/motion data is missing';
          else
            run(i).info.complete='yes';
            run(i).info.validwindow=[0 run(i).data_nirs.data_raw.time{1}(end)];
          end
        end
    end