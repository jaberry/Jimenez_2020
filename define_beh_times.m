
function [stim_times_matrix,numb_stim,stim_lengths]= define_beh_times(type,FPS)
%this function will generate a cell array with the row indices for
%different session behavior stimuli. These indices will be
%used in other functions for running analyses in diff behavioral states
%specified
%'CFCD1SHOCK' for the 2-second shock in CFCD1 encoding session
%'CFCD2SHOCK' for the (3) 2-second tone-shocks after the CFCD2 retreival session
%'CFCD2TONES' for the (3) 18-second tones (pre tone-shocks) after the CFCD2 retreival session
%'CFCD3TONES' for the (3) 20-second tones after the CFCD3 novel Context B session

if strcmpi(type,'CFCD1SHOCK')==1
    numb_stim=2;
    stim1=1:180*FPS; %D1 first 180s
    stim2=(180*FPS)+1:182*FPS;%D1 shock
    stim_times_matrix{numb_stim}=NaN;
    
    stim_times_matrix={stim1, stim2};
    stim_lengths(numb_stim)=NaN;
    stim_lengths=horzcat(length(stim1), length(stim2));
    
elseif strcmpi(type,'CFCD2SHOCK')==1
    numb_stim=2;
    stim1=1:180*FPS; %D2 first 180s
    stim2=[198*FPS+1:200*FPS,258*FPS+1:260*FPS,318*FPS+1:320*FPS];%D2 all 3 tone-shocks
    stim_times_matrix{numb_stim}=NaN;
    
    stim_times_matrix={stim1, stim2};
    stim_lengths(numb_stim)=NaN;
    stim_lengths=horzcat(length(stim1), length(stim2));
    
    elseif strcmpi(type,'CFCD2TONES')==1
    numb_stim=2;
    stim1=1:180*FPS; %D2 first 180s
    stim2=[(180*FPS)+1:198*FPS,(240*FPS)+1:258*FPS,(300*FPS)+1:318*FPS];%D2 all 3 tones (pre-shock)
    stim_times_matrix{numb_stim}=NaN;
    
    stim_times_matrix={stim1, stim2};
    stim_lengths(numb_stim)=NaN;
    stim_lengths=horzcat(length(stim1), length(stim2));
    
    elseif strcmpi(type,'CFCD3TONES')==1
    numb_stim=2;
    stim1=1:180*FPS; %D3 first 180s
    stim2=[(180*FPS)+1:200*FPS,(240*FPS)+1:260*FPS,(300*FPS)+1:320*FPS];%D3 all 3 tones, full time
    stim_times_matrix{numb_stim}=NaN;
    
    stim_times_matrix={stim1, stim2};
    stim_lengths(numb_stim)=NaN;
    stim_lengths=horzcat(length(stim1), length(stim2));
    
end