function [zscored_cell,cell_transients,cell_events, cell_AUC]=detect_ca_transients(raw_cell,thresh,baseline,t_half,FR)
%Ca transient event detection, adapted from Dombeck 2007
%author Jessica Jimenez, Columbia University, jcj2123@columbia.edu

%inputs
%raw_cell= Ca2 transient data in TxN format, N=cells (columns), T=time (rows), raw format
%thresh= minimum amplitude size of ca transient in s.d.
%baseline= s.d. offset value of ca transient
%t_half= half life of gcamp type used (s), if taken from Chen TW et al Nature 2013, gcamp6f t_half=0.200 s
%FR= framerate of raw_cell

%outputs
%cell_transients= TxN matrix with s.d. values for all the timepoints of the qualified transients (all zeros except for transients, red trace in fig)
%cell_events= TxN matrix with calcium transient peak events (all zeros except for the amplitude value at the peak timepoint of each transient, asterisk in fig)
%cell_AUC= TxN matrix with calcium transient area under the curve (AOC) values (all zeros except for the AOC value assigned to the peak timepoint of each transient)
%zscored_cell= TxN matrix with zscored raw_cell data (blue trace in fig)

% the above 4 outputs will be saved into a .mat in your directory called
% "ca2events", and each cell figure will also be saved showing the event
% detection performance

%zscore your data
pophist=reshape(raw_cell,[],1); % generate single vector with all cell fluorescence values
pop_offset=quantile(pophist,0.50); %find the 50% quantile value (for "silent" time points)
silent=pophist<pop_offset; %find timepoints without ca transients based on threshold above
mu=mean(pophist(silent==1)); % specify mu from the silent timepoints
[~, ~, sigma] = zscore(raw_cell,1); %specify sigma from the entire time series
zscored_cell = bsxfun(@rdivide, bsxfun(@minus, raw_cell, mu), sigma); %convert transients into zscores using mu from silent timepoints and sigma from all timepoints

celldata=zscored_cell;

%preallocate outputs
tk=size(celldata);
cell_transients=zeros(tk);
cell_events=zeros(tk);
cell_AUC=zeros(tk);

%define minimum duration of calcium transient based on gcamp type used
decayrate=0.693/t_half; %simplified from (-ln(A/Ao)/t_half), [A/Ao]=0.5 at t half-life, [-ln(A/Ao)]=0.693
minduration=-(log(baseline/thresh))/decayrate; %minimum (s) duration for ca transient of minimum specified s.d. amplitude threshold
minframes= round(minduration*FR); %minimum number of frames the ca transient should last above baseline

%identify qualified ca transients and generate outputs
for k=1:size(celldata,2);
    
    onset=find(celldata(:,k)>thresh); %find all timepoints where flourescence greater than threshold
    offset=find(celldata(:,k)>baseline); %find all timepoints where floursecence greater than baseline (transient offset)
    
    found=1;
    for m = 1:length(offset)-1
        
        if found == 1
            start = offset(m); %specify start index of transient from offset vector
            found = 0;
        end
        if offset(m+1) - offset(m) > 1 %specify stop index of transient from offset vector
            finish = offset(m);
            [M,I]=max(celldata(start:finish,k)); %find the peak value in that start-stop range
            transientvect=start:finish;
            maxamp_ind=transientvect(I); %retrieve "cell" index of that peak value
            peak_to_offset_vect=maxamp_ind:finish;
            found  = 1;
            
            if ismember(maxamp_ind,onset)>0 && length(peak_to_offset_vect)>=minframes; %if the peak value index from start-stop in offset is also found in onset vector, the transient exceeded the 2SD threshold
                cell_transients(start:finish,k) = celldata(start:finish,k); %retrieve "cell" values for all the timepoints of that transient
                cell_events(maxamp_ind,k)=M; %create a matrix with all the calcium transient peak events (all zeros except for the amplitude value at the peak timepoint)
                transient_area=trapz(celldata(start:finish,k)); %integrate the area under the curve of the transient from start-stop
                cell_AUC(maxamp_ind,k)=transient_area; %create a matrix with all the calcium transient AOC values (all zeros except for the AOC value assigned to the peak timepoint)
            end
            
        end
        if m== length(offset)-1 %dealing with the last index in the vector, same as above
            finish= offset(m+1);
            [M,I]=max(celldata(start:finish,k));
            transientvect=start:finish;
            maxamp_ind=transientvect(I);
            peak_to_offset_vect=maxamp_ind:finish;
            found  = 1;
            if ismember(maxamp_ind,onset)>0  && length(peak_to_offset_vect)>=minframes;
                cell_transients(start:finish,k) = celldata(start:finish,k);
                cell_events(maxamp_ind,k)=M;
                transient_area=trapz(celldata(start:finish,k));
                cell_AUC(maxamp_ind,k)=transient_area;
            end
        end
        
    end
    
end




% detect multi-peak transients & update in cell_events

for k=1:size(cell_transients,2);
    
    [~, time]=findpeaks(cell_transients(:,k),'MinPeakProminence',1.5,'MinPeakDistance',FR); %built-in matlab 'findpeaks' fxn
    %minpeak distance is 1 sec (as specified by your frame rate), and min peak
    %prominence must be 1.5SD in size
    cell_events(time,k)=cell_transients(time,k); % cell events with multi-peak transients added
end


%plot and save figures with event detection results
X=(1:size(cell_events,1))';
events=cell_events;
events(events==0)=NaN;

for i=1:size(cell_events,2);
    figure(i);plot(zscored_cell(:,i),'b'); hold on; plot(cell_transients(:,i),'r'); hold on; plot(X,events(:,i),'m*');title(['Cell',num2str(i)],'Fontsize',10);
end

h = get(0,'children');
length_h=1:length(h);
length_h=sort(length_h,'descend');
for i=1:length(h)
    saveas(h(length_h(i)), ['figure' num2str(i)], 'fig');
end

save('ca2events.mat','cell_transients', 'cell_events', 'cell_AUC', 'zscored_cell')

end

