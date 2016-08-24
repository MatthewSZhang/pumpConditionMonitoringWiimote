% Copied from XL spreadsheet
conditionVec = [4
6
2
6
7
7
2
3
8
6
1
1
9
8
9
9
1
7
1
5
1
10
1
9
7
3
7
2
1
2]';
[counts,centers] = hist(conditionVec);
figure;bar(centers,counts);
grid on;
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = {'1 Working excellent',...
'2 Noisy but working',...
'3 Dry borehole',...
'4 Rising main leak',...
'5 Broken Seal',...
'6 Unsure',...
'7 Water leaking from pump',...
'8 Worn bush bearing',...
'9 Stiff handle',...
'10 Worn seal'};
title('Distribution of "condition" of pumps','FontSize',14);