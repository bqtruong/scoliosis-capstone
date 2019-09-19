load('full_stats.mat');
ages_str = a;
ages = zeros(length(ages_str),1);
for x = 1:length(ages_str)
    if ages_str{x}(end) == 'M'
        ages(x) = str2double(ages_str{x}(1:3)) / 12;
    elseif ages_str{x}(end) == 'Y'
        ages(x) = str2double(ages_str{x}(1:3));
    else
        disp('Something is wrong!');
    end
end 
a = ages;

% Convert mm^3 to L
v = v .* 1e-6;

% Outlier adjustment
a116 = a == 116;
h(a116) = [];
w(a116) = [];
v(a116) = [];
g(a116) = [];
a(a116) = [];

hw0 = h == 0 | w == 0;
v(hw0) = [];
g(hw0) = [];
a(hw0) = [];
h(hw0) = [];
w(hw0) = [];

% Messed up MS
v(36) = [];
g(36) = [];
a(36) = [];
h(36) = [];
w(36) = [];

g = cellfun(@(x) x == 'M', g);
in = [a, g, h, w]';
out = v';
corrs = zeros(1,100);
k = 10;
indices = crossvalind('Kfold', out, 10);

rmses = zeros(1,k);
maes = zeros(1,k);
maes2 = zeros(1,k);
for i = 1:k
    in_train_cross = in(:, indices == i);
    out_train_cross = out(indices == i);
    net = fitnet(3,'trainbr');
    net = train(net, in_train_cross, out_train_cross);
    in_test_cross = in(:, indices ~= i);
    out_test_cross = out(indices ~= i);
    out_test_net = net(in_test_cross);
    rmses(i) = sqrt(mean((out_test_cross - out_test_net).^2));
    maes(i) = mean(abs(out_test_cross - out_test_net));
    maes2(i) = mae(out_test_cross - out_test_net);
end

% rmses_all = [];
% rmses_e_all = [];
% maes_all = [];
% maes_e_all = [];
% for neuron_count = 2:10
%     rmses = zeros(1,k);
%     maes = zeros(1,k);
%     maes2 = zeros(1,k);
%     for i = 1:k
%         in_train_cross = in(:, indices == i);
%         out_train_cross = out(indices == i);
%         net = fitnet(neuron_count,'trainbr');
%         net = train(net, in_train_cross, out_train_cross);
%         in_test_cross = in(:, indices ~= i);
%         out_test_cross = out(indices ~= i);
%         out_test_net = net(in_test_cross);
%         rmses(i) = sqrt(mean((out_test_cross - out_test_net).^2));
%         maes(i) = mean(abs(out_test_cross - out_test_net));
%         maes2(i) = mae(out_test_cross - out_test_net);
%     end
%     disp(['MAE ' num2str(neuron_count) ' ' num2str(mean(maes)) ' ' num2str(std(maes))]);
%     disp(['RMSE ' num2str(neuron_count) ' ' num2str(mean(rmses)) ' ' num2str(std(rmses))]);
%     maes_all = [maes_all mean(maes)];
%     maes_e_all = [maes_e_all std(maes)];
%     rmses_all = [rmses_all mean(rmses)];
%     rmses_e_all = [rmses_e_all std(rmses)];
% end

% net = fitnet(3,'trainbr');
% net = train(net,in,out);
% y = net(in);
% corr = corrcoef(out,y);
% corrs = corr(1,2) .^ 2;
% e = out-y;
% MAE = mae(e);
% disp(mean(corrs));
% plotregression(out, y, 'All');
% h = gcf;
% a=findobj(h,'type','axe');
% xl=get(get(a,'xlabel'),'string');
% yl=get(get(a,'ylabel'),'string');
% yl = strrep(yl, 'Output', 'Predicted');
% yl = strrep(yl, 'Target', 'Actual');
% yl = strrep(yl, '~', '');
% xlabel('Actual Volume (L)');
% ylabel('Predicted Volume (L)');
% title('Linear regression of the network');
% text(.5,.1,['r^{2} = ' num2str(corr(1,2).^2)],'Units','normalized');
% text(.5,.05,yl,'Units','normalized');

% maes = zeros(1,100);
% stderror = zeros(1,100);
% parfor i = 1:100
%     net = fitnet(3,'trainbr');
%     net = train(net,in,out);
%     y = net(in);
%     corr = corrcoef(out,y);
%     e = out-y;
%     stderror(i) = sqrt(sum(e .^ 2) / length(e));
%     MAE = mae(e);
% end
% disp(stderror);

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
