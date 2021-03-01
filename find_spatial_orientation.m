function find_spatial_orientation(ellipsoid)

[U, ~, ~]   = svd(ellipsoid);
U  = U(:, 1:10); 

U(abs(U) > 1)  = sign(U(abs(U) > 1));

axis_angle  = acos(U) * 180 / pi;
axis_angle(axis_angle > 90)  = 180 - axis_angle(axis_angle > 90);

figure, imagesc(axis_angle, [0, 90]), colorbar, title('Orientation')
figure, plot(min(axis_angle)), title('Orientation')

end