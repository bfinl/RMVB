function my_ellipsoid(h, A, O, n)
[X, Y, Z]  = sphere(n); points  = cat(3, X, Y, Z);
for i = 1 : n+1
    for j = 1 : n+1
        points(i, j, :) = (A * squeeze(points(i, j, :)) + O).';
    end
end
surf(points(:, :, 1), points(:, :, 2), points(:, :, 3), 'facealpha', 0.2) 
end