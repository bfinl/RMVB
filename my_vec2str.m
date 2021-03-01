function str = my_vec2str(vec)

str = [];
for i = 1 : length(vec)
    if(i == 1), str = [str, sprintf('%.2f', vec(i))];
    else, str = [str, sprintf(',\t%.2f', vec(i))]; end
end

end