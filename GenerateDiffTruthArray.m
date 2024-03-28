%% Generate lin fits with 95% confidence interval and truth array
%
% Written by Mark I. Grimes and Mick D. Mantle, University of Cambridge, March 2024.
%
% Data should be inputted as an (n x [m*2]) matrix called "all_data", where n is the number of concentrations, and m is the number of stressed fractions studied. 
% Data should be put in positions (:,1:m); errors should be put in positions (:,m:end), following same order as data.
%
% Aggregate content values should be inputted as an (n x m) matrix called "AggVals", where n is the number of concentrations, and m is the number of stressed
% fractions studied.
% 

critical_value = 1.96; % change as req'd; chosen for 95% confidence interval generation

xRows = 1;
disp(' ');
xColumns = input('Please input the number of concentrations you have: '); % User to input concentrations studied
disp(' ');

disp('Input concentration values: ');
Concs = zeros(xRows,xColumns);
for k = 1:xRows
    for m = 1:xColumns
        Concs(k,m) = input("Input the matrix value for (" + k + "," + m +"): ");
    end
end

f1 = figure;
count1 = 0;
y_pred=[];x_pred=[];
Markers = {'ok','dk','sk','^k','pk'}; % for marker style - change as req'd
MarkerFC = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880]}; % for marker colours - change as req'd

for j = 1:size(all_data,1)
    
x = AggVals(j,:);

% Scatter plot
DataPlots(j) = plot(x, all_data(1 + count1, 1:6), Markers{j}, 'MarkerSize', 12,'LineWidth',1.5,'MarkerFaceColor',MarkerFC{j});
 
hold on;
 
% Error bars
errorbar(x, all_data(1 + count1, 1:6), all_data(1 + count1, 7:end), 'k', 'LineStyle', 'None', 'CapSize', 5,'linewidth',1.5);
 
% Linear regression
coefficients = polyfit(x, all_data(1 + count1, 1:6), 1);
line = polyval(coefficients, x);
plot(x, line, 'Color', 'b', 'LineWidth',1.5);
 
% Confidence intervals
x_pred = x;
y_pred(j,:) = polyval(coefficients, x_pred);
 
result = mean(all_data(:,7:end), 'all');
result_array = zeros(6) + result;

y_err = result_array;  % Assuming the error bars represent the standard error of the mean
y_pred_upper = y_pred(j,:) + critical_value * y_err(j); % Generate values for confidence interval fill area
y_pred_lower = y_pred(j,:) - critical_value * y_err(j);
 
% Fill between confidence intervals
fill([x_pred, fliplr(x_pred)], [y_pred_lower, fliplr(y_pred_upper)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');


count1 = count1 + 1;
end
 
hold off;
 
% Plot data
legend([DataPlots],num2str(Concs','%5.2f mg mL^{-1}'),'Location','EastOutside'); % change unit as req'd
set(gca,'FontSize',20,'FontWeight','bold','LineWidth',2,'box','on');
xlabel('% aggregate','FontSize',30)
ylabel('{\it D}(^1H_2O) [m^2 s^{-1}]', 'FontSize',30)
f1.Position = [400 50 1200 750];
set(gcf, 'renderer', 'painters')
ylim([(min(min(all_data(:,1:6)))*0.997) (max(max(all_data(:,1:6)))*1.003)]) % set as req'd
xLimBuffer = max(max(AggVals))*1.05 - max(max(AggVals));
xlim([(min(min(AggVals)))-xLimBuffer max(max(AggVals))+xLimBuffer])
yticks([2.57e-9 2.58e-9 2.59e-9 2.6e-9 2.61e-9 2.62e-9]) % set as req'd
axis square
ytickformat('%.2f')

%% Set up +/- confidence intervals array will be 2m x n in size
% Also set up +/- D values array of same dimensions to see when values are equal

D = input("Input the diffusion value of interest: ");  % User to input diffusion coefficient value of interest
D_array = zeros(1,6) + D;  % make a 1D array of 6 identical values

result = mean(all_data(:, 7:end), 'all');  % calculates average of errors
result_array = zeros(1,6) + result;  % repeats single value 6 times

D_upper = D_array + result_array;  % upper value of D array
D_lower = D_array - result_array;  % lower value of D array

y_pred_L = [];  % create dummy cell array to store upper 95% conf
y_pred_U = [];  % create dummy cell array to store lower 95% conf

y_err2 = result_array;  % set y_error to result array (is a single value)
ci_array = zeros(10, 6);
d_array = zeros(10, 6);
count1 = 1;
y_pred_upper2 = y_pred + (critical_value * y_err2);
y_pred_lower2 = y_pred - (critical_value * y_err2);

for j = 1:2:(size(all_data,1)*2)
    y_pred_U = y_pred_upper2(count1,:);  % stores this upper value for one BSA conc
    ci_array(j,:) = y_pred_U;
    d_array(j,:) = D_upper;

    count1 = count1 + 1;
end

count1 = 1;
for k = 2:2:(size(all_data,1)*2)
    y_pred_L = y_pred_lower2(count1,:);  % stores this Lower value for one BSA conc
    ci_array(k,:) = y_pred_L;
    d_array(k,:) = D_lower;
    count1 = count1 + 1;
end

% Set the threshold for comparison (0.1%)
threshold = 0.001;

% Perform element-wise comparison
truth_array = abs(ci_array - d_array) <= threshold*abs(ci_array);

% Display the comparison array in command window
disp('Comparison Array:');
disp(truth_array);

% Plot the comparison array
f2 = figure;
y_ticks = 1.5:2:9.5;
x_ticks = 0.5:1:5.5;
y_labels = Concs;
x_labels = [];
imagesc(truth_array);
hold on
yline(2.5,'Color','w','LineWidth',5,'Alpha',1); yline(4.5,'Color','w','LineWidth',5,'Alpha',1) % lines for formatting; change number if req'd
yline(6.5,'Color','w','LineWidth',5,'Alpha',1); yline(8.5,'Color','w','LineWidth',5,'Alpha',1)
title("Truth array for {\it D} = " + (D*10^9) + " \times 10^{-9} m^2 s^{-1}");
set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels, 'YTick', y_ticks,'YTickLabel',y_labels);
set(gca,'FontSize',20,'FontWeight','bold','LineWidth',2,'box','on');
ylabel('Concentration [mg mL^{-1}]', 'FontSize',30)
f2.Position = [400 50 950 750];
set(gcf, 'renderer', 'painters')
axis square
colormap(parula(2));
clb = colorbar;
clb.Ticks = [0.25 0.75];
clb.TickLabels = [0 1];


