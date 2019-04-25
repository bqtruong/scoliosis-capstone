load('full_stats.mat');
ages_str = a;
ages = zeros(length(ages_str),1);
for x = 1:length(ages_str)
    if ages_str{x}(end) == 'M'
        ages(x) = str2double(ages_str{x}(1:3)) / 12;
    elseif ages_str{x}(end) == 'Y'
        ages(x) = str2double(ages_str{x}(1:3));
    else
        disp("Something's wrong!");
    end
end 
a = ages;
% Outlier adjustment
h(a == 116) = [];
w(a == 116) = [];
v(a == 116) = [];
g(a == 116) = [];
a(a == 116) = [];
h(40) = [];
w(40) = [];
v(40) = [];
g(40) = [];
a(40) = [];
g = cellfun(@(x) x == 'M', g);
in = [a, g, h, w]';
out = v';
net = fitnet(3,'trainbr');
net = train(net,in,out);
y = net(in);
corr = corrcoef(out,y);
e = out-y;
MAE = mae(e);

%% Testing best number of neurons
% neurons = 2:25;
% ntrials = 10;
% all_corrs = zeros(ntrials, length(neurons));
% parfor trials = 1:ntrials
%     corrs = zeros(1,length(neurons));
%     for i = 1:length(neurons)
%         net = fitnet(i,'trainbr');
%         net = train(net,in,out);
%         y = net(in);
%         corr = corrcoef(out, y);
%         corrs(i) = corr(1,2);
%     end
%     all_corrs(trials, :) = corrs;
% end
% all_corrs_mean = mean(all_corrs, 1);
% plot(neurons, all_corrs_mean);
