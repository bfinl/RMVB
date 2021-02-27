function h = error_bar(Y, STAT)

Color(1, :)  = [0.1289, 0.4844, 0.6602];
Color(2, :)  = [0.4258, 0.7344, 0.8828];
Color(3, :)  = [0.9570, 0.4609, 0.0664];
Color(4, :)  = [1.0000, 0.8500, 0.0000];

Y_mid  = zeros(size(Y));
U_err  = zeros(size(Y));
L_err  = zeros(size(Y));

for i = 1 : size(Y, 1)
    for j = 1 : size(Y, 2)
        pro_Y  = Y{i, j};
        pro_Y  = pro_Y(~isnan(pro_Y));
%         if( isempty(pro_Y) )
%             pro_Y  = -1;
%         end
        
        switch STAT
            case 'mean'
                Y_mid(i, j)  = mean(pro_Y);
                U_err(i, j)  = +std(pro_Y, 1)/2;
                L_err(i, j)  = -std(pro_Y, 1)/2;  
                
            case 'median'
                Y_mid(i, j)  = median(pro_Y);
                U_err(i, j)  = quantile(pro_Y, 0.75)-Y_mid(i, j);
                L_err(i, j)  = Y_mid(i, j)-quantile(pro_Y, 0.25);
        end
    end
end

numbars = size(Y, 2);
groupwidth = min(0.8, numbars/(numbars+1.5));
X  = 1 * (1:size(Y, 1));
h  = bar(X, Y_mid, 1, 'FaceColor', 'flat');

count = 1;
for i = 1 : length(h)
    set(h(i),'FaceColor', Color(count, :))
    count = count + size(Color, 1)/length(h);
end

% set(h(1),'FaceColor', Color1)
% set(h(2),'FaceColor', Color2)
% set(h(3),'FaceColor', Color3)
% set(h(4),'FaceColor', Color4)

% set(h(1),'FaceColor', Color1)
% set(h(2),'FaceColor', Color3)

% set(h(1),'FaceColor', Color1)

hold on
for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = X - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, Y_mid(:, i), L_err(:, i), U_err(:, i), 'k', 'linestyle', 'none');
end


end
