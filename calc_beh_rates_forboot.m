
function[eventrates_bybeh,states,stim_length]=calc_beh_rates_forboot(session,type,FPS)
%outputs rates and events per cell per behavior state specified
%calls on the 'define_beh_times' function to define behavior time stamps for
%diff imaging sessions
%inputs: 
%session=cell data, and specify what imaging session as
%type= specify the behavior session as one of the below:
%'CFCD1SHOCK' for the 2-second shock in CFCD1 encoding session
%'CFCD2SHOCK' for the (3) 2-second tone-shocks after the CFCD2 retreival session
%'CFCD2TONES' for the (3) 18-second tones (pre tone-shocks) after the CFCD2 retreival session
%'CFCD3TONES' for the (3) 20-second tones after the CFCD3 novel Context B session

[states,numb_stim,stim_length]= define_beh_times(type,FPS);
%uses 'define_beh_times' function to define behavior columns to analyze
%outputs: the behavior indices in a cell array "states", how many behavior
%states are extracted from that imaging session (numb_beh), and the length of time the
%animal spent in that behavior state to use for rate calculations
%(beh_length)
%eventrates_bybeh: column 1= rates in cond1, column 2= rates in cond2,
%column 3= rates cond2-rates cond1, column4= rates cond1-rates cond2

k=1:size(session,2);
eventrates_bybeh(k,4)=zeros;
cellevents_behstates(k,4)=zeros;
for b=1:2
    events_bycell=sum(session(states{b},:),1);
    %loops through each behavioral state in cell array to sum events per cell in that state
    cellevents_behstates(k,b)=events_bycell';
    %cells in rows, # of events in each behavior, each behavior in diff column
end

for b=1:numb_stim
    eventrates_bybeh(k,b)=cellevents_behstates(:,b)./stim_length(b);
    %calculates rate by dividing number of events by the length of that
    %behavioral state
end
eventrates_bybeh(k,3)=eventrates_bybeh(k,2)-eventrates_bybeh(k,1);
%cond 2 selective: rate differences column2-column1, to be used in bootstrap script later

eventrates_bybeh(k,4)=eventrates_bybeh(k,1)-eventrates_bybeh(k,2);
%cond1selective: rate differences column1-column2, to be used in bootstrap script later

end