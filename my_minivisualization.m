function [intersectionmap, ratiomap, anglemap] = my_minivisualization(LFD, A, d, disp_opt, disp_dim, disp_space, varargin)
%%
[Uproj, ~, ~]  = svds(disp_space, disp_dim); 
display('=================================================================%')
me  = (norm(LFD-Uproj*Uproj.'*LFD, 'fro')/norm(LFD, 'fro'))^2;
fprintf('The missing energy for the purpose of 3D display: %%%f\n', me)
display('=================================================================%')
h = figure; grid, hold on, axis equal, xlabel('x'), ylabel('y'), zlabel('z')


%%
for i = 1 : d : size(LFD, 2)
    ctr  = Uproj.' * LFD(:, i:i+d-1);
    scatter3(ctr(1, :), ctr(2, :), ctr(3, :), 1000 *ones(1, d), '*', 'LineWidth', 2)
    if(i == 1), scatter3(ctr(1, :), ctr(2, :), ctr(3, :), 1000 *ones(1, d), 'o', 'LineWidth', 2), end
    for r = 1 : d    
        switch disp_opt
            case 'projection' % The projection of the high dimensional space onto the low dimensional display space
                P  = A(:, :, i+r-1) * (A(:, :, i+r-1).'); 
                [U, S, ~]  = svds(Uproj * (Uproj.') * P * Uproj * (Uproj.'), disp_dim);
                Aprospace  = Uproj.' * U * (S^0.5);
                
            case 'intersection' % The intersection of the high dimensional and low dimensional display spaces
                iP  = (A(:, :, i+r-1) * (A(:, :, i+r-1).'))^-1; 
                [U, S, ~]  = svds(Uproj * (Uproj.') * iP * Uproj * (Uproj.'), disp_dim);  
                Aprospace  = Uproj.' * U * (S^-0.5); 
        end
        my_ellipsoid(h, Aprospace , ctr(:, r), 20)
        mArrow3(ctr(:, r).', ctr(:, mod(r, d)+1).', 'color', 'red', 'stemWidth', 0.001, 'facealpha', 1);
        mArrow3(zeros(1, d), ctr(:, r).', 'color', 'black', 'stemWidth', 0.005, 'facealpha', 1);       
    end
end
for i = 1 : size(LFD, 2)/d
    if(~isempty(varargin))
      INDINF   = varargin{2};
      LFD_NBR  = varargin{1};
      20*log10(norm(LFD_NBR(:, :, i), 'fro')/norm(Uproj*Uproj.'*LFD_NBR(:, :, i) - LFD_NBR(:, :, i), 'fro'))
      scatter3(Uproj(:, 1).'*LFD_NBR(:, :, i), Uproj(:, 2).'*LFD_NBR(:, :, i), Uproj(:, 3).'*LFD_NBR(:, :, i), 10*ones(1, size(LFD_NBR, 2)), 'ro', 'LineWidth', 2)
      scatter3(Uproj(:, 1).'*LFD_NBR(:, INDINF{i}, i), Uproj(:, 2).'*LFD_NBR(:, INDINF{i}, i), Uproj(:, 3).'*LFD_NBR(:, INDINF{i}, i), 100*ones(1, size(INDINF{i}, 2)), 'k*', 'LineWidth', 2)
    end
end


%%
for i = 1 : d : size(LFD, 2)
    for r = 1 : d
        [U, S, ~]   = svd(A(:, :, i+r-1));
        for s = 1 : d
            dir  = LFD(:, i+s-1) - LFD(:, i+r-1);
            dis  = norms(dir); dir  = dir / dis;
            rad1 = find_ellipradius(A(:, :, i+r-1), dir);
            rad2 = find_ellipradius(A(:, :, i+s-1), dir);
            strongestaxis_angle  = acos(U(:, 1).' * dir) * 180 / pi;
            strongestaxis_angle(strongestaxis_angle > 90)  = 180 - strongestaxis_angle(strongestaxis_angle > 90);
            intersectionmap(r, s, i) = (dis - rad1 - rad2) / dis;
            anglemap(r, s, i)  = strongestaxis_angle;
            ratiomap(r, s, i)  = rad1 / S(1);       
        end
    end
end


end








% function [intersectionmap, ratiomap, anglemap] = my_minivisualization(lfd, A, d, disp_type, varargin)
% 
% for r = 1 : d
%     [U, S, ~]   = svd(A(:, :, r));
%     for s = 1 : d
%         dir  = lfd(:, s) - lfd(:, r);
%         dis  = norms(dir); dir  = dir / dis;
%         rad1 = find_ellipradius(A(:, :, r), dir);
%         rad2 = find_ellipradius(A(:, :, s), dir);
%         strongestaxis_angle  = acos(U(:, 1).' * dir) * 180 / pi;
%         strongestaxis_angle(strongestaxis_angle > 90)  = 180 - strongestaxis_angle(strongestaxis_angle > 90);
%         intersectionmap(r, s) = (dis - rad1 - rad2) / dis;
%         anglemap(r, s)  = strongestaxis_angle;
%         ratiomap(r, s)  = rad1 / S(1);       
%     end
% end
% 
% if(strcmp(disp_type, 'Display'))
%     [Upro, ~, ~]  = svds(lfd, d); ctr  = Upro.' * lfd;
%     h = figure; scatter3(ctr(1, :), ctr(2, :), ctr(3, :), 1000 *ones(1, d), '*', 'LineWidth', 2)
%     hold on, axis equal, xlabel('x'), ylabel('y'), zlabel('z')
%     for r = 1 : d        
%         iP  = (A(:, :, r) * (A(:, :, r).'))^-1; 
%         [U, S, ~]  = svds(Upro * (Upro.') * iP * Upro * (Upro.'), d);   
%         Arot(:, :, r) = Upro.' * U * (S^-0.5);
%         my_ellipsoid(h, Arot(:, :, r) , ctr(:, r), 20)
%         mArrow3(ctr(:, r).', ctr(:, mod(r, d)+1).', 'color', 'red', 'stemWidth', 0.001, 'facealpha', 1);
%         mArrow3(zeros(1, d), ctr(:, r).', 'color', 'black', 'stemWidth', 0.001, 'facealpha', 1);       
%     end
%     if(~isempty(varargin))
%       INDINF   = varargin{2};
%       LFD_NBR  = varargin{1};
%       20*log10(norm(LFD_NBR, 'fro')/norm(Upro*Upro.'*LFD_NBR - LFD_NBR, 'fro'))
%       scatter3(Upro(:, 1).'*LFD_NBR, Upro(:, 2).'*LFD_NBR, Upro(:, 3).'*LFD_NBR, 10*ones(1, size(LFD_NBR, 2)), 'ro', 'LineWidth', 2)
%       scatter3(Upro(:, 1).'*LFD_NBR(:, INDINF), Upro(:, 2).'*LFD_NBR(:, INDINF), Upro(:, 3).'*LFD_NBR(:, INDINF), 100*ones(1, size(INDINF, 2)), 'k*', 'LineWidth', 2)
%     end
% end
% 
% end