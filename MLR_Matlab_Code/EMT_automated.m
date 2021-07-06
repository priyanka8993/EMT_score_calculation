%GPL570

load('RelevantData.mat');

dirInfo = importdata("matlab_input_temp.txt")
 
gseData = importdata(char(dirInfo(1,1)), '\t', 1);
data = gseData.data(1:3,:);
norm_Val = gseData.data(6:25,:);

%Pre-normalization scatterplots
%1. %(CLDN7, VIM/CDH1)
figure; hold on;
plot(0,0,'.k'); plot(0,0,'.b');
scatter(DataNCI60EMTGenes8(b1,:),DataNCI60EMTGenes8(b2,:),'k*');
b3=1; b4=2; b5=3;
scatter(data(b3,:),data(b4,:)./data(b5,:),'b*');
xlabel('CLDN7 Expression'); ylabel('VIM/CDH1 Expression');
legend('NCI60','DataGSE');

nci_indices = importdata(char(dirInfo(2,1)));
nci_indices = nci_indices(6:25);
NCI_data = DataNCI60(nci_indices,:);
NCIMean = nanmean(mean(NCI_data));

mean_val= mean(mean(norm_Val));
d=mean_val-NCIMean; 
NormDataGSE=data-d;




%Post-normalization scatterplots
%1. (CLDN7, VIM/CDH1)figure; hold on;
plot(0,0,'.k'); plot(0,0,'.b');
scatter(DataNCI60EMTGenes8(b1,:),DataNCI60EMTGenes8(b2,:),'k*');
scatter(data(b3,:),data(b4,:)./data(b5,:),'b*');
scatter(NormDataGSE(b3,:),NormDataGSE(b4,:)./NormDataGSE(b5,:),'r*');
xlabel('CLDN7 Expression'); ylabel('VIM/CDH1 Expression');
legend('NCI60','data','NormDataGSE');


%saveas(gcf,fullfile(char(dirInfo(1,1)),'png');


%Predictions
%1. (CLDN7, VIM/CDH1)
[YhatDataGSEEMTGenes8pair1, PredictionsDataGSEEMTGenes8pair1] = MLR3(...
    B1, GeneList1,...
    [{'CLDN7'}; {'VIM/CDH1'}],...
    [data(b3,:); data(b4,:)./data(b5,:)]);

[YhatDataGSEEMTGenes8pair1Norm, PredictionsDataGSEEMTGenes8pair1Norm] = MLR3(...
    B1, GeneList1,...
    [{'CLDN7'}; {'VIM/CDH1'}],...
    [NormDataGSE(b3,:); NormDataGSE(b4,:)./NormDataGSE(b5,:)]);

NUM = size(data, 2);
ScoreEMT3 = nan(NUM, 1);
for j = 1:NUM
    if YhatDataGSEEMTGenes8pair1(j,1) > YhatDataGSEEMTGenes8pair1(j,3)
        ScoreEMT3(j) = YhatDataGSEEMTGenes8pair1(j, 2);
    else
        ScoreEMT3(j) = 2.0 - YhatDataGSEEMTGenes8pair1(j, 2);
    end
end

ScoreEMT3_norm = nan(NUM,1);
for j = 1:NUM
    if YhatDataGSEEMTGenes8pair1Norm(j, 1) > YhatDataGSEEMTGenes8pair1Norm(j, 3)
        ScoreEMT3_norm(j) = YhatDataGSEEMTGenes8pair1Norm(j, 2);
    else
        ScoreEMT3_norm(j) = 2.0 - YhatDataGSEEMTGenes8pair1Norm(j, 2);
    end
end

labels = gseData.textdata(1,2:end);
sampleName = transpose(labels);
T = table(sampleName,ScoreEMT3_norm);
writetable(T,char(dirInfo(3,1)),'Delimiter','\t','WriteRowNames',false);