function setup_display(CDR, currytri, LOC, DPYLOC, TH, ttl, varargin)

    n = size(LOC, 2);
    m = size(DPYLOC, 2); 
    nTim = size(CDR, ndims(CDR));

    V  = zeros(m+n, nTim); 
    if(ndims(CDR) == 2)
        V(1:n, 1:nTim) = CDR;
    elseif(ndims(CDR) == 3)
        V(1:n, 1:nTim) = squeeze(CDR(4, :, :));
    end
    V  = abs(V);

    if(isempty(varargin))
        epl = {}; 
    else
        epl  = varargin{1};
        epl_disp_type = varargin{2};
        if(~iscell(epl)), epl = {epl}; end
    end

    curryloc  = [LOC DPYLOC];

    
    %%
    nMin = min(curryloc,[],2)-10;
    nMax = max(curryloc,[],2)+10;
    h_fig  = figure;
    h_fig  = ApplyProperties(h_fig, 'Customized-v1'); 
    handles.ptr = 1;
    guidata(h_fig, handles);
    ttl(ttl == '_')  = '-';
    axis([nMin(1), nMax(1), nMin(2), nMax(2), nMin(3), nMax(3)]);
    axis equal; axis vis3d;
%     view(0, 0);
%     view( 107, 29);
    view(215, 45);
%     cmp  = colormap(colorcube); 
    cmp  = colormap(hot); 
    cmp  = cmp(1:48, :);
    colormap(cmp)
    colorbar
    
    
    %%
    rsctn_display()
    set(h_fig, 'KeyPressFcn', @(h_obj, evt) KeyPress(h_obj, evt));
    mini_display(V(:, 1), {ttl, ['time: ' num2str(1)], 'Mode: Custom'})

    function KeyPress(h_obj, evt)

    %     handles  = guidata(h_obj);

        if(strcmp(evt.Key, 'rightarrow') || strcmp(evt.Key, 'uparrow'))
            handles.ptr  = handles.ptr + 1;
            handles.ptr  = handles.ptr - (handles.ptr > nTim);
            mini_display(V(:, handles.ptr), {ttl, ['time: ' num2str(handles.ptr)], 'Mode: Custom'})

        elseif(strcmp(evt.Key, 'leftarrow') || strcmp(evt.Key, 'downarrow'))
            handles.ptr  = handles.ptr - 1;
            handles.ptr  = handles.ptr + (handles.ptr < 1);
            mini_display(V(:, handles.ptr), {ttl, ['time: ' num2str(handles.ptr)], 'Mode: Custom'})

        elseif(strcmp(evt.Key, 'p') && ismember('shift', get(gcbo,'currentModifier')))
            mini_display(max(V, [], 2), {ttl, ['time: ', 'NA'], 'Mode: Peak'})
            
        elseif(strcmp(evt.Key, 'w') && ismember('shift', get(gcbo,'currentModifier')))
            mini_display(norms(V, [], 2), {ttl, ['time: ', 'NA'], 'Mode: Power'})
            
        elseif(strcmp(evt.Key, 's') && ismember('shift', get(gcbo,'currentModifier')))
            mini_display(std(V, 1, 2), {ttl, ['time: ', 'NA'], 'Mode: STD'})

        elseif(strcmp(evt.Key, 'a') && ismember('shift', get(gcbo,'currentModifier')))
            mini_display(mean(V, 2), {ttl, ['time: ', 'NA'], 'Average Mode'}) 

        elseif(strcmp(evt.Key, 'm') && ismember('shift', get(gcbo,'currentModifier')))
            for ii = 1 : nTim
                mini_display(V(:, ii), {ttl, ['time: ', num2str(ii)], 'Mode: Movie'}) 
                pause(0)
            end
            handles.ptr  = nTim;    
        elseif(strcmp(evt.Key, 'r') && ismember('shift', get(gcbo,'currentModifier')))
            mini_display([], 'Mode: Resection') 
            rsctn_display()
        end

        guidata(h_obj, handles);

    end

    function mini_display(V, ttl)
        title(ttl)
        mn = -eps; mx = eps; 
        mp = 0.8906*ones(size(curryloc, 2), 3); 
        if(~isempty(V)), [mp, mn, mx]  = find_map(V, mp); end
        hpatch = patch('vertices', curryloc', 'faces', currytri', 'FaceVertexCData', mp);
        set(hpatch, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceLighting', 'phong', 'SpecularStrength', 0.2,...
        'DiffuseStrength', 0.8, 'AmbientStrength', 0.3, 'FaceAlpha', 1);
        delete(findall(gcf,'Type','light'))
        camlight right;
        caxis([mn, mx])         
    end

    function [mp, mn, mx] = find_map(V, mp)
        mx  = max(V);
        ig  = find(V >= mx*TH);
%         ig  = 1:length(V);
        mn  = min(V(ig));
        if(mx == mn), mx  = mn + mn*1e-6; end
        mp(ig, :)  = cmp(floor((length(cmp)-1)*(V(ig)-mn)/(mx-mn)+1), :); 
    end

    function rsctn_display()
        hold on
        for i_epl = 1 : length(epl)
            rsctnloc  = epl{i_epl};
            if (size(rsctnloc, 2) == 1)
                scatter3(rsctnloc(1), rsctnloc(2), rsctnloc(3), 500, 'co', 'LineWidth', 2);
            else
                switch epl_disp_type
                    case 'pnt' % for point sources
                        scatter3(rsctnloc(1, :), rsctnloc(2, :), rsctnloc(3, :), 500, 'co', 'LineWidth', 2);
                    case 'dot'
                        scatter3(rsctnloc(1, :), rsctnloc(2, :), rsctnloc(3, :), 1, 'c.', 'LineWidth', 2);
                    case 'net'
                        [~, rsctnind] = find_nvoxel(rsctnloc, curryloc); rsctnind = unique(rsctnind);
                        rsctnloc  = curryloc(:, rsctnind);
                        rsctntri  = currytri(:, find_triind(currytri, rsctnind, 3));
                        for ii = 1 : size(rsctnloc, 2)
                            rsctntri(rsctntri == rsctnind(ii))  = ii;
                        end
                        patch('vertices', rsctnloc', 'faces', rsctntri', 'EdgeColor', [0 100 220]/256, 'FaceColor', [0 100 220]/256);
                        delete(findall(gcf,'Type','light'))
                        camlight right;
                    case 'boundary'
                        [~, rsctnind] = find_nvoxel(rsctnloc, curryloc); rsctnind = unique(rsctnind);
                        rsctnloc  = curryloc(:, rsctnind);
                        rsctntri  = currytri(:, find_triind(currytri, rsctnind, 3));
                        for ii = 1 : size(rsctnloc, 2)
                            rsctntri(rsctntri == rsctnind(ii))  = ii;
                        end
                        bnd_edge  = find_boundary_edges(rsctntri);
                        line([rsctnloc(1, bnd_edge(1, :)); rsctnloc(1, bnd_edge(2, :))],...
                             [rsctnloc(2, bnd_edge(1, :)); rsctnloc(2, bnd_edge(2, :))],...
                             [rsctnloc(3, bnd_edge(1, :)); rsctnloc(3, bnd_edge(2, :))], 'LineWidth', 2.5, 'Color', [0 100 220]/256 ...
                            )
                        
                        
                    otherwise
                        disp('Wrong resection display mode!!!')
                end
            end
        end
    end


end

