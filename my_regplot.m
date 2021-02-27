function my_regplot(costx, costy, Alpha, Lambda)
                
h = figure; 
ApplyProperties(h, 'Customized-v1');

hold on
for i = 1 : size(costx, 1)
%     subplot(1, size(cost1, 1), i), 
    plot(costx(i, :), costy(i, :));  axis equal
    lgnd{i}  = ['\alpha = ' num2str(Alpha(i))];
end
hold off
xlabel('Data Fitting Term'), ylabel('Regularization')
title(['\alpha = ' num2str(Alpha(i))]), legend(lgnd);
datatipalpha  = repmat(Alpha, [1, length(Lambda)]);
datatiplambda = repmat(Lambda, [length(Alpha), 1]);
datatiplambda = datatiplambda(:).'; 
dcm  = datacursormode(h); datacursormode on
set(dcm, 'updatefcn', {@myregparamdatatip, costx(:), costy(:), datatipalpha, datatiplambda})
                
end

function output_txt = myregparamdatatip(obj, event_obj, X, Y, Alpha, Lambda)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
output_txt = {['X: ', num2str(pos(1),4)],...
              ['Y: ', num2str(pos(2),4)]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end

ind  = find(X==pos(1) & Y==pos(2));          
output_txt{end+1} = ['\Alpha: ', num2str(Alpha(ind),4)];
output_txt{end+1} = ['\Lambda: ', num2str(Lambda(ind),4)];

end