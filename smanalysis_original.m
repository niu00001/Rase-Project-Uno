%% Stock Analysis Algorithm

% Analyzes stock using various methods including NASDAQ Dozen. 
% Based on data extracted from finviz. 

% Ideas:
% Include volume into the analysis. Overbought/oversold.
% Use of next earning date info.

% SEC website to extract more data

clear
clc

%% Data Sheet Input & Processing

rawdata = readtable('SP500_20171202.xlsx','ReadVariableNames',false);
datafile =  table2cell(rawdata);

ndatafile = datafile;
ndata_L = length(ndatafile);

for k = 2:ndata_L
   for j = [7:68,70:size(datafile,2)]

       incomp=datafile{k,j};
       incomp=strrep(incomp,'B','e9');
       incomp=strrep(incomp,'M','e6');
       ndatafile{k,j}=incomp;

   end
end

data_num = str2double(ndatafile);

%% Information Indexing

% Revenue Info (1)

GrossM_ind = find(strcmp(ndatafile(1,:),'Gross M')==1);
ProfitM_ind = find(strcmp(ndatafile(1,:),'Profit M')==1);

PerfQ_ind = find(strcmp(ndatafile(1,:),'Perf Quart')==1);
PerfHY_ind = find(strcmp(ndatafile(1,:),'Perf Half')==1);
PerfY_ind = find(strcmp(ndatafile(1,:),'Perf Year')==1);
PerfYTD_ind =  find(strcmp(ndatafile(1,:),'Perf YTD')==1);

data_num(data_num(:,PerfQ_ind)>0,PerfQ_ind)=1;
data_num(data_num(:,PerfQ_ind)<0,PerfQ_ind)=0;

data_num(data_num(:,PerfHY_ind)>0,PerfHY_ind)=1;
data_num(data_num(:,PerfHY_ind)<0,PerfHY_ind)=0;

data_num(data_num(:,PerfY_ind)>0,PerfY_ind)=1;
data_num(data_num(:,PerfY_ind)<0,PerfY_ind)=0;

data_num(data_num(:,PerfYTD_ind)>0,PerfYTD_ind)=1;
data_num(data_num(:,PerfYTD_ind)<0,PerfYTD_ind)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Earnings Per Share Info (2)

EPSQtoQ_ind = find(strcmp(ndatafile(1,:),'EPS Q/Q')==1);
EPSthisY_ind = find(strcmp(ndatafile(1,:),'EPS this Y')==1);
EPSpast5Y_ind = find(strcmp(ndatafile(1,:),'EPS past 5Y')==1);

data_num(data_num(:,EPSQtoQ_ind)>0,EPSQtoQ_ind)=1;
data_num(data_num(:,EPSQtoQ_ind)<0,EPSQtoQ_ind)=0;

data_num(data_num(:,EPSthisY_ind)>0,EPSthisY_ind)=1;
data_num(data_num(:,EPSthisY_ind)<0,EPSthisY_ind)=0;

data_num(data_num(:,EPSpast5Y_ind)>0,EPSpast5Y_ind)=1;
data_num(data_num(:,EPSpast5Y_ind)<0,EPSpast5Y_ind)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Return of Equity (3)

ROE_ind = find(strcmp(ndatafile(1,:),'ROE')==1);

data_num(data_num(:,ROE_ind)>0,ROE_ind)=1;
data_num(data_num(:,ROE_ind)<0,ROE_ind)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analyst Recommendations (4)

Recom_ind = find(strcmp(ndatafile(1,:),'Recom')==1);

data_num(data_num(:,Recom_ind)<=3,Recom_ind)=1;
data_num(data_num(:,Recom_ind)>3,Recom_ind)=0;

% finviz uses following scale:
% 1 = strong buy
% 2 = buy
% 3 = hold
% 4 = sell
% 5 = strong sell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Earnings Surprise (5)

% finviz has no info on this. Check if you can extract from Nasdaq. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Earnings Forecast (6)

EPSnextY_ind = find(strcmp(ndatafile(1,:),'EPS next Y')==1);

data_num(data_num(:,EPSnextY_ind)>0,EPSnextY_ind)=1;
data_num(data_num(:,EPSnextY_ind)<0,EPSnextY_ind)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Earnings Growth (7)

EPSnext5Y_ind = find(strcmp(ndatafile(1,:),'EPS next 5Y')==1);

data_num(data_num(:,EPSnext5Y_ind)>0,EPSnext5Y_ind)=1;
data_num(data_num(:,EPSnext5Y_ind)<0,EPSnext5Y_ind)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PEG Ratio (8)

PEG_ind = find(strcmp(ndatafile(1,:),'PEG')==1);

data_num(data_num(:,PEG_ind)<1,PEG_ind)=1;
data_num(data_num(:,PEG_ind)>1,PEG_ind)=0;

% PEG less than 1 is favorable, PEG greater than 1 is not favorable.
% Realistically its more of a scale since its calculated with the following
% formula PEG = (P/E)/Annual EPS Growth. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Industry Price-Earnings (9)

PE_ind = find(strcmp(ndatafile(1,:),'P/E') == 1);
Sector_ind = find(strcmp(ndatafile(1,:),'Sector') == 1);
Tech_ind = find(strcmp(ndatafile(:,Sector_ind),'Technology') == 1);
Serv_ind = find(strcmp(ndatafile(:,Sector_ind),'Services') == 1);
BM_ind = find(strcmp(ndatafile(:,Sector_ind),'Basic Materials') == 1);
Fin_ind = find(strcmp(ndatafile(:,Sector_ind),'Financial') == 1);
CG_ind = find(strcmp(ndatafile(:,Sector_ind),'Consumer Goods') == 1);
Health_ind = find(strcmp(ndatafile(:,Sector_ind),'Healthcare') == 1);
IG_ind = find(strcmp(ndatafile(:,Sector_ind),'Industrial Goods') == 1);
Uti_ind = find(strcmp(ndatafile(:,Sector_ind),'Utilities') == 1);

TechPEAvg = nanmean(data_num(Tech_ind, PE_ind));
TechPEMed = nanmedian(data_num(Tech_ind, PE_ind));
TechPEStd = nanstd(data_num(Tech_ind, PE_ind));

ServPEAvg = nanmean(data_num(Serv_ind, PE_ind));
ServPEMed = nanmedian(data_num(Serv_ind, PE_ind));
ServPEStd = nanstd(data_num(Serv_ind, PE_ind));

BMPEAvg = nanmean(data_num(BM_ind, PE_ind));
BMPEMed = nanmedian(data_num(BM_ind, PE_ind));
BMPEStd = nanstd(data_num(BM_ind, PE_ind));

FinPEAvg = nanmean(data_num(Fin_ind, PE_ind));
FinPEMed = nanmedian(data_num(Fin_ind, PE_ind));
FinPEStd = nanstd(data_num(Fin_ind, PE_ind));

CGPEAvg = nanmean(data_num(CG_ind, PE_ind));
CGPEMed = nanmedian(data_num(CG_ind, PE_ind));
CGPEStd = nanstd(data_num(CG_ind, PE_ind));

HealthPEAvg = nanmean(data_num(Health_ind, PE_ind));
HealthPEMed = nanmedian(data_num(Health_ind, PE_ind));
HealthPEStd = nanstd(data_num(Health_ind, PE_ind));

IGPEAvg = nanmean(data_num(IG_ind, PE_ind));
IGPEMed = nanmedian(data_num(IG_ind, PE_ind));
IGPEStd = nanstd(data_num(IG_ind, PE_ind));

UtiPEAvg = nanmean(data_num(Uti_ind, PE_ind));
UtiPEMed = nanmedian(data_num(Uti_ind, PE_ind));
UtiPEStd = nanstd(data_num(Uti_ind, PE_ind));



% % compare to average
% data_num(data_num(Tech_ind,PE_ind) < TechPEAvg, PE_ind) = 1;
% data_num(data_num(Serv_ind,PE_ind) < ServPEAvg, PE_ind) = 1;
% data_num(data_num(BM_ind,PE_ind) < BMPEAvg, PE_ind) = 1;
% data_num(data_num(Fin_ind,PE_ind) < FinPEAvg, PE_ind) = 1;
% data_num(data_num(CG_ind,PE_ind) < CGPEAvg, PE_ind) = 1;
% data_num(data_num(Health_ind,PE_ind) < HealthPEAvg, PE_ind) = 1;
% data_num(data_num(IG_ind,PE_ind) < IGPEAvg, PE_ind) = 1;
% data_num(data_num(Uti_ind,PE_ind) < UtiPEAvg, PE_ind) = 1;
% 
% data_num(data_num(Tech_ind,PE_ind) > TechPEAvg, PE_ind) = 0;
% data_num(data_num(Serv_ind,PE_ind) > ServPEAvg, PE_ind) = 0;
% data_num(data_num(BM_ind,PE_ind) > BMPEAvg, PE_ind) = 0;
% data_num(data_num(Fin_ind,PE_ind) > FinPEAvg, PE_ind) = 0;
% data_num(data_num(CG_ind,PE_ind) > CGPEAvg, PE_ind) = 0;
% data_num(data_num(Health_ind,PE_ind) > HealthPEAvg, PE_ind) = 0;
% data_num(data_num(IG_ind,PE_ind) > IGPEAvg, PE_ind) = 0;
% data_num(data_num(Uti_ind,PE_ind) > UtiPEAvg, PE_ind) = 0;

%compare to median
data_num(data_num(Tech_ind,PE_ind) < TechPEMed, PE_ind) = 1;
data_num(data_num(Serv_ind,PE_ind) < ServPEMed, PE_ind) = 1;
data_num(data_num(BM_ind,PE_ind) < BMPEMed, PE_ind) = 1;
data_num(data_num(Fin_ind,PE_ind) < FinPEMed, PE_ind) = 1;
data_num(data_num(CG_ind,PE_ind) < CGPEMed, PE_ind) = 1;
data_num(data_num(Health_ind,PE_ind) < HealthPEMed, PE_ind) = 1;
data_num(data_num(IG_ind,PE_ind) < IGPEMed, PE_ind) = 1;
data_num(data_num(Uti_ind,PE_ind) < UtiPEMed, PE_ind) = 1;

data_num(data_num(Tech_ind,PE_ind) > TechPEMed, PE_ind) = 0;
data_num(data_num(Serv_ind,PE_ind) > ServPEMed, PE_ind) = 0;
data_num(data_num(BM_ind,PE_ind) > BMPEMed, PE_ind) = 0;
data_num(data_num(Fin_ind,PE_ind) > FinPEMed, PE_ind) = 0;
data_num(data_num(CG_ind,PE_ind) > CGPEMed, PE_ind) = 0;
data_num(data_num(Health_ind,PE_ind) > HealthPEMed, PE_ind) = 0;
data_num(data_num(IG_ind,PE_ind) > IGPEMed, PE_ind) = 0;
data_num(data_num(Uti_ind,PE_ind) > UtiPEMed, PE_ind) = 0;

PE_per_sector=[TechPEMed, ServPEMed, BMPEMed, FinPEMed, CGPEMed, HealthPEMed, IGPEMed, UtiPEMed];
Sector_names={ 'Tech'; 'Serv'; 'Basic Materials'; 'Financial'; 'Consumer Goods'; 'Health'; 'Industrial Goods'; 'Utilities'};
disp(['The sector with lowest median PE is' Sector_names(PE_per_sector==min(PE_per_sector))])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Days to Cover (10)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Insider Trading (11)

InsiderTrans_ind = find(strcmp(ndatafile(1,:),'Insider Trans')==1);

data_num(data_num(:,InsiderTrans_ind)>0,InsiderTrans_ind)=1;
data_num(data_num(:,InsiderTrans_ind)<0,InsiderTrans_ind)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sector Analysis (12)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Relative Strength Index (13)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Forward PE (14)

FPE_ind = find(strcmp(ndatafile(1,:),'Fwd P/E') == 1);

data_num(data_num(:,FPE_ind) < data_num(:,PE_ind),FPE_ind) = 1;
data_num(data_num(:,FPE_ind) > data_num(:,PE_ind),FPE_ind) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Stock Score

ind=[PerfQ_ind,PerfHY_ind,PerfY_ind,PerfYTD_ind,...
    EPSQtoQ_ind,EPSthisY_ind,EPSnext5Y_ind,ROE_ind,Recom_ind,...
    EPSnextY_ind,EPSnext5Y_ind,PEG_ind,InsiderTrans_ind, PE_ind,FPE_ind];

stock_aggregate = data_num(:,ind);
   
%% Stock Analyzer

stock_score = nansum(stock_aggregate,2);

maxval_gold = max(stock_score)
idx_gold = find(stock_score == maxval_gold);
idx_silver = find(stock_score == (maxval_gold-1));

gold = string(ndatafile(idx_gold,3))
silver = string(ndatafile(idx_silver,3))