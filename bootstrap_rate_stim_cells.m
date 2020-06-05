
function [boot_rates, signif_thresholds, signif_cell_inds,signif_cell_totals]= bootstrap_rate_stim_cells (session, type, threshold,FPS, iter)
%this fxn identifies selective cells based on the difference in rate
%between 2 behavioral states (ie open-closed arms EPM/EZM, center-periphery
%OFT, novelob-oppfield novel object task etc)
%calls on "shake" function from matlab file exchange https://www.mathworks.com/matlabcentral/fileexchange/10067-shake
%INPUTS:
%session= ca event data matrix, rows=time, columns=cells
%type= specify the behavior session as one of the below:
%'CFCD1SHOCK' for the 2-second shock in CFCD1 encoding session
%'CFCD2SHOCK' for the (3) 2-second tone-shocks after the CFCD2 retreival session
%'CFCD2TONES' for the (3) 18-second tones (pre tone-shocks) after the CFCD2 retreival session
%'CFCD3TONES' for the (3) 20-second tones after the CFCD3 novel Context B session
%FPS= frames per second
%signif threshold= specify significance of cell's true ca2+ rate from bootrap distribution
%'all'= 1SD & 2SD cut offs, '1SD', or '2SD' only
%iter= number of iterations for bootstrap

%OUTPUTS:
%boot_rates= 1x2 cell array with the rates generated from the bootstrap, element 1 is condition1 boot rates, element 2 for condition 2. Within each element is a matrix of rates organized by rows=cells, columns= iterations
%signif thresholds= ca2+ event rate values at significant thresholds for
%each cell
%signif_cell_inds= cell indices of selective cells (non-selective, cond1
%selective, then cond2 selective) for 'all', a 1x6 cell array (first for
%1SD, then for 2SD), otherwise a 1x3 cell array
%signif_cell_totals= vector with the total number of selective cells, same
%format as cell array for signif_cell_inds

    session=session>0; %make sure the ca event data matrix is a logical (amplitude values removed)
    
    [cellrates_behstates,states,beh_length]=calc_beh_rates_forboot(session,type,FPS);
    %pre-allocate bootstrap output sizes
    k=1:size(session,2);
    cond1_events(k,iter)=zeros;
    cond2_events(k,iter)=zeros;
for i=1:1000%number of shuffle iterations
             cellshuffle=shake(session,1); %make sure the 'shake' fxn is in your matlab path, shuffles the cell event data
             cond1_events_temp = cellshuffle(states{1},:); 
             %makes a temp matrix w the shuffled cell data w only the condition 1 rows
             cond1_events(:,i)=sum(cond1_events_temp);
             %sums the # of events per cell from that temp matrix
             
             %same as above but for condition2
            cond2_events_temp = cellshuffle(states{2},:);
             cond2_events(:,i)=sum(cond2_events_temp);
end

cond1_rate=cond1_events./beh_length(1); %converts events to rate
cond2_rate=cond2_events./beh_length(2); 
   
    boot_rates={cond1_rate,cond2_rate};
    %outputs a cell array w all of the bootstrap rates from above for future
    %reference
    
    signif_thresh_cond2_2SD= quantile(cond2_rate,0.95,2); %2SD quantile from bootstrap rate
    signif_thresh_cond2_1SD= quantile(cond2_rate,0.68,2); %1SD " "
    
    signif_thresh_cond1_2SD= quantile(cond1_rate,0.95,2);
    signif_thresh_cond1_1SD= quantile(cond1_rate,0.68,2);
    
    signif_thresholds=horzcat(signif_thresh_cond1_1SD, signif_thresh_cond2_1SD, signif_thresh_cond1_2SD,signif_thresh_cond2_2SD);
    %outputs the signif threshold levels per cell for future reference
    
    cond2_selective_cells_2SD=cellrates_behstates(:,2)>signif_thresh_cond2_2SD;
    %logical for if cell's diff rate (cond2) exceeds signif threshold
    cond2_selective_cells_1SD=cellrates_behstates(:,2)>signif_thresh_cond2_1SD;
    
    cond1_selective_cells_2SD=cellrates_behstates(:,1)>signif_thresh_cond1_2SD;
    %logical for if cell's diff rate (cond1) exceeds signif threshold
    cond1_selective_cells_1SD=cellrates_behstates(:,1)>signif_thresh_cond1_1SD;

selectivity_logicalindex=horzcat(cond1_selective_cells_1SD,cond2_selective_cells_1SD,cond1_selective_cells_2SD,cond2_selective_cells_2SD);

%%%gets cell indices in each condition
cond2_selective_cells2SD=find(cond2_selective_cells_2SD>0);
cond2_selective_cells1SD=find(cond2_selective_cells_1SD>0);

cond1_selective_cells2SD=find(cond1_selective_cells_2SD>0);
cond1_selective_cells1SD=find(cond1_selective_cells_1SD>0);

nonselective_cells2SD=find(cond2_selective_cells_2SD<1 & cond1_selective_cells_2SD<1);
nonselective_cells1SD=find(cond2_selective_cells_1SD<1 & cond1_selective_cells_1SD<1);
%%%

selective_cell_indices_1SD={nonselective_cells1SD, cond1_selective_cells1SD, cond2_selective_cells1SD};
selective_cell_indices_2SD={nonselective_cells2SD, cond1_selective_cells2SD, cond2_selective_cells2SD};

selective_cell_totals_1SD=horzcat(length(nonselective_cells1SD), length(cond1_selective_cells1SD), length(cond2_selective_cells1SD));
selective_cell_totals_2SD=horzcat(length(nonselective_cells2SD), length(cond1_selective_cells2SD), length(cond2_selective_cells2SD));

if strcmpi(threshold,'all')==1
    signif_cell_inds=horzcat(selective_cell_indices_1SD,selective_cell_indices_2SD);
    signif_cell_totals=horzcat(selective_cell_totals_1SD,selective_cell_totals_2SD);
elseif strcmpi(type,'1SD')==1
    signif_cell_inds=selective_cell_indices_1SD;
    signif_cell_totals=selective_cell_totals_1SD;
elseif strcmpi(type,'2SD')==1
    signif_cell_inds=selective_cell_indices_2SD;
    signif_cell_totals=selective_cell_totals_2SD;
end
save('rate_selective_cells.mat','signif_cell_totals','signif_cell_inds','signif_thresholds','boot_rates');
end