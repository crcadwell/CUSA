cd /Users/cathryn/Dropbox/'CUSA DP Project'/Code
addpath venn

%% Figure 2 panels C - D

clear all; close all
load data.mat

figure;

X = [length(find(score==1)) length(find(score==2)) length(find(score==3)) length(find(score==4))];
% X = [113 23 108 52]
labels = {'Complete scan ','Beyond coverslip','Missing tissue <0.5 mm','Missing tissue >=0.5 mm'};
subplot(1,3,1)
pie(X,labels);

d = find(score==3);
%venn([sum(ROI(d)) sum(MSA(d))],length(find(ROI(d)+MSA(d)==2)))
% Values are 105, 10, 7

d = find(score==4);
%venn([sum(ROI(d)) sum(MSA(d))],length(find(ROI(d)+MSA(d)==2)))
% Values are 45, 13, 6

d = find(score>=3);
subplot(1,3,2)
venn([sum(ROI(d)) sum(MSA(d))],length(find(ROI(d)+MSA(d)==2)))
% Values are 150, 23, 13  
ax = gca;
ax.XLim = [-8 10];
ax.Visible = 'off';
axis equal
legend({'ROI','MSA'})
legend('boxoff')

% Comparison of MSA contribution between WSI with small v. large missing
% tissue
X = [1 2];
[Y(1),ci(1,:)] = binofit(3, 108, 0.1);
[Y(2),ci(2,:)] = binofit(7, 52, 0.1);
subplot (1,3,3)
scatter(X,Y,'k'); hold on;
errorbar(X,Y,Y'-ci(:,1),ci(:,2)-Y','k');
ax = gca;
ax.XLim = [0 3];
ax.Children(1).LineStyle = 'none'
[h,p,stats] = fishertest([105 3; 45 7])
%p=0.0141

%Comparison of ROI contribution between WSI with small v. large missing
%tissue
[h,p,stats] = fishertest([97 11; 39 13])
%p = 0.0185

% Comare rates of infidelity for 2016 v. 2018 slides
t1 = find(year==2016); t2 = find(year==2018);
d1 = find(score>=3);
d2 = find(score==3);
d3 = find(score==4);

X1 = [length(intersect(t1,d1)) length(t1)-length(intersect(t1,d1));...
    length(intersect(t2,d1)) length(t2)-length(intersect(t2,d1))];
[h,p,stats] = fishertest(X1)
%p=0.4853 for all infidelities

X2 = [length(intersect(t1,d2)) length(t1)-length(intersect(t1,d2));...
    length(intersect(t2,d2)) length(t2)-length(intersect(t2,d2))];
[h,p,stats] = fishertest(X2)
%p=0.0906 for small infidelities

X3 = [length(intersect(t1,d3)) length(t1)-length(intersect(t1,d3));...
    length(intersect(t2,d3)) length(t2)-length(intersect(t2,d3))];
[h,p,stats] = fishertest(X3)
%p=0.2223 for large infidelities


%% Figure 3 panels

clear all; close all
load data.mat

figure;

for i = 1:length(score)
    if strcmp(ID_on_WSIS{i},'Yes') && score(i)>=3
        if strcmp(ID_on_glass{i},'Yes') && score(i)>=3
            ID_score(i) = 4;
        else
            ID_score(i) = 2;
        end
    elseif strcmp(ID_on_glass{i},'Yes') && score(i)>=3
        ID_score(i) = 3;
    else
        ID_score(i) = 1;
    end
end

X = [length(find(ID_score==1)) length(find(ID_score==2)) length(find(ID_score==4)) length(find(ID_score==3))];
% X = [136 66 63 31]
labels = {'Complete scan ','Infidelity on WSIS only','Infidelity on WSIS, additional on glass','Infidelity on glass only'};
subplot(1,3,1)
pie(X,labels);
ax1 = gca;

X = [1 2 3 5 6 7];
s = find(score==3);
d(1,1) = length(intersect(s, find(ID_score==2)));
d(2,1) = length(intersect(s, find(ID_score==4)));
d(3,1) = length(intersect(s, find(ID_score==3)));
l = find(score==4);
d(4,1) = length(intersect(l, find(ID_score==2)));
d(5,1) = length(intersect(l, find(ID_score==4)));
d(6,1) = length(intersect(l, find(ID_score==3)));
d(1:3,2) = length(s);
d(4:6,2) = length(l);
for i = 1:size(d,1)
    [Y(i),ci(i,:)] = binofit(d(i,1), d(i,2), 0.1);
end
subplot (1,3,2:3)
for i=1:3
    scatter(X([i i+3]),Y([i i+3]),'filled'); hold on;
    errorbar(X([i i+3]),Y([i i+3]),Y([i i+3])'-ci([i i+3],1),ci([i i+3],2)-Y([i i+3])','k');
end
ax = gca;
ax.XLim = [0 8];
ax.Children(1).LineStyle = 'none'; 
ax.Children(3).LineStyle = 'none'; 
ax.Children(5).LineStyle = 'none'; 
legend(ax.Children([6 4 2]), {'Infidelity on WSIS only','Infidelity on WSIS, additional on glass',...
    'Infidelity on glass only'});
legend('boxoff')

X2table = [d(1:3,1)';d(4:6,1)'];
[chi2stat,pOverall] = ChiSquared(X2table)
%chi2stat = 7.5742
%p = 0.0227
%df = (3-1)*(2-1) = 2

[~,pPostHoc(1)] = ChiSquared([X2table(:,1) sum(X2table(:,2:3),2)]);
[~,pPostHoc(2)] = ChiSquared([X2table(:,2) sum(X2table(:,[1 3]),2)]);
[~,pPostHoc(3)] = ChiSquared([X2table(:,3) sum(X2table(:,1:2),2)]);
pPostHoc = pPostHoc*3 %Bonferroni correction
%pPostHoc = [0.1711 2.5682 0.0284];


