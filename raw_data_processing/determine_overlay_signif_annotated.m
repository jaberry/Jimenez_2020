%% determine if the overlap between 2 selective cell populations is greater/less than you would expect by chance

%% create a mock distribution with following characteristics:
%	- vector of zeros with length (m) where m= true number of cells analyzed across mice
%   - '1's randomly assigned to (y) elements where y= # of
%       selective cells in condition 1

m=848; %total cells analyzed
y=111; %total number of condition 1 selective cells 
mock_distribution=zeros(m,1); %create a vector of zeros, size of total # cells analyzed
selective_cells=randperm(m,y); %randomly select indices for y number of cells from the total
mock_distribution(selective_cells,1)=1; %assign those randomly selected cells a value of 1 to represent "selective cell"
assign_check=sum(mock_distribution); %should be equal to y
%% randomly select (n) # of cells from the mock distribution, 10,000 times

n=78; %number of selective cells in condition 2
overlapcells(10000,1)=NaN; %pre-allocate
for i=1:10000
    
    selectivecells=randsample(mock_distribution,n);
    overlapcells(i,1)=sum(selectivecells);
end

%convert # overlap cell rand sample output to "% of condition 2"
mock_overlap_percentofcond2=(overlapcells./n).*100;
%% determine pvalue of true overlap

mean=mean(mock_overlap_percentofcond2);
sigma=std(mock_overlap_percentofcond2);
true_overlap_percentofcond2=(34/78)*100; %true % overlap % of cond1/cond2 
[h,p,ci,zval]=ztest(true_overlap_percentofcond2,mean,sigma);

statsconcat=horzcat(true_overlap_percentofcond2,h,p,mean-(sigma*2),mean+(sigma*2),zval,ci(1),ci(2));
