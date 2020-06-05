function [corrpairs,corrpairs_byFOV,D1_corr, D2_corr, D3_corr,D1_corr_ind,D2_corr_ind,D3_corr_ind,D1_corr_values,D2_corr_values,D3_corr_values]=correlated_pairs_Rthresh(m,session,FPS,session_length,time_bin,Rthresh)
%calculates numb of correlated pairs per cell given 0.3 threshold as per
%Vogel nature paper

%INPUTS:
%m= mouse number
%session= a 1x3 cell array with Ca2+ event data matrices (rows=time, columns= cells); {1}=CFC day1, {2}=CFC day2, {3}=CFC day3 
%FPS= frames per second
%session_length= time in seconds
%time_bin= length of time over which to bin the Ca2+ event data (in
%seconds)
%Rthresh= pearson's R threshold used to define a correlated pair (R=0.3 for
%the paper)

%OUTPUTS (all saved into matlab path under "mouse#_corr_activity")
%corrpairs= matrix with numb correlated pairs/ cell (rows=cells, columns=
%session days)
%corrpairs_byFOV= numb corrpairs per cell/ total # cells in FOV
%D1_corr, D2_corr, and D3_corr= matrix for each session day with pearson's R values for each
%cell combination (cellxcell matrix of R values)

%D1_corr_ind, D2_corr_ind, D3_corr_ind= across 3 days, cell array with vectors of cell
%indices for each cell's significantly correlated cell partners (rows= corrpair indices for each
%cell). The first value in each vector is the cellID for that cell. 
%For example, if cell 3 is correlated with cell 18, 52, and 68 on CFC day2, then in D2_corr_ind, at the 3rd row, the cell array
%will contain the following vector: [3,18,52,68]

%D1_corr_values, D2_corr_values, D3_corr_values= accross 3 days, cell array with vectors of R values for each cell's significantly correlated cell pair (rows= Rvalues for each
%cell). 
%continuing example of 'cell 3' above, on day 2, if Rvalue=0.41 for cell 3 to cell 18, 0.32 for cell 3 to cell 52, and 0.67 for cell 3 to cell 68, then D2_corr_values will contain the following vector at the 3rd row: [0.41, 0.32, 0.67]


    N=size(session{1},2); %number of cells
 
        session1=session{1};
        session1=session1(1:session_length*FPS,:)>0; 
        session2=session{2};
        session2=session2(1:session_length*FPS,:)>0;
        session3=session{3};
        session3=session3(1:session_length*FPS,:)>0;
        
        frames_perbin=FPS*time_bin; %how many frames will be summed per bin
        b=size(session1,1)/frames_perbin;%number of bins in time window
        
        %bin the data into the time window desired (ie 1 sec, 3 sec etc)
        k=1:size(session1,2);
        time_window_sum1(b,k)=zeros;
        for k=1:size(session1,2)
            time_window1=reshape(session1(:,k),frames_perbin,b);
            time_window_sum1(:,k)=sum(time_window1,1); %#events per cell over those binned frames
        end
        
        time_window_sum2(b,k)=zeros;
        for k=1:size(session2,2)
            time_window2=reshape(session2(:,k),frames_perbin,b);
            time_window_sum2(:,k)=sum(time_window2,1);
        end
        
        time_window_sum3(b,k)=zeros;
        for k=1:size(session3,2)
            time_window3=reshape(session3(:,k),frames_perbin,b);
            time_window_sum3(:,k)=sum(time_window3,1);
        end
        
        D1_corr=corrcoef(time_window_sum1); %correlation between every cell pair in that binned cell data
        D2_corr=corrcoef(time_window_sum2);
        D3_corr=corrcoef(time_window_sum3);
        
        D1_corr_logical=(D1_corr>Rthresh & D1_corr<1); %which cell pairs are correlated above the 0.3 significance level
        D1_numbpairs=sum(D1_corr_logical,2); %number of correlated pairs per cell
        for i=1:size(D1_corr_logical,2)
            D1_corr_ind{i,1}=i; %define index cell
            D1_corr_ind{i,1}=vertcat(D1_corr_ind{i,1},find(D1_corr_logical(:,i)==1)); %cell ID of index cell, followed by cell ID of each of its correlated pairs
            D1_corr_values{i,1}=D1_corr(D1_corr_logical(:,i)==1,i); %pearson's R values for each correlated pair per cell
        end

        D2_corr_logical=(D2_corr>Rthresh & D2_corr<1);
        D2_numbpairs=sum(D2_corr_logical,2);
        for i=1:size(D2_corr_logical,2)
            D2_corr_ind{i,1}=i;
            D2_corr_ind{i,1}=vertcat(D2_corr_ind{i,1},find(D2_corr_logical(:,i)==1));
            D2_corr_values{i,1}=D2_corr(D2_corr_logical(:,i)==1,i);
        end
        
        D3_corr_logical=(D3_corr>Rthresh & D3_corr<1);
        D3_numbpairs=sum(D3_corr_logical,2);
        for i=1:size(D3_corr_logical,2)
            D3_corr_ind{i,1}=i;
            D3_corr_ind{i,1}=vertcat(D3_corr_ind{i,1},find(D3_corr_logical(:,i)==1));
            D3_corr_values{i,1}=D3_corr(D3_corr_logical(:,i)==1,i);
        end
      
        corrpairs(k,3)=zeros;
        corrpairs=horzcat(D1_numbpairs,D2_numbpairs,D3_numbpairs); %corr pairs per cells across days
        corrpairs_byFOV=corrpairs./N;

    file_name = sprintf('mouse%d_corr_activity.mat',m);
    save(file_name,'corrpairs','corrpairs_byFOV','D1_corr','D2_corr','D3_corr','D1_corr_ind','D2_corr_ind','D3_corr_ind','D1_corr_values','D2_corr_values','D3_corr_values');

end
