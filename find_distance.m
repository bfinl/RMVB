function dis  = find_distance(src_loc, dst_loc, opt)

switch opt
    
    case 'Max-Min'
        for i_src = 1 : size(src_loc, 2)
            dis(i_src) = min(norms(repmat(src_loc(:, i_src), [1, size(dst_loc, 2)]) - dst_loc));            
        end
        dis  = max(dis);
        
        
    case 'Ave-Min'
        for i_src = 1 : size(src_loc, 2)
            dis(i_src) = min(norms(repmat(src_loc(:, i_src), [1, size(dst_loc, 2)]) - dst_loc));            
        end
        dis  = mean(dis);
        
        
    case 'Max-Ave'
        for i_src = 1 : size(src_loc, 2)
            dis(i_src) = mean(norms(repmat(src_loc(:, i_src), [1, size(dst_loc, 2)]) - dst_loc));            
        end
        dis  = max(dis);
        
        
    case 'Min'
        for i_src = 1 : size(src_loc, 2)
            dis(i_src) = min(norms(repmat(src_loc(:, i_src), [1, size(dst_loc, 2)]) - dst_loc));            
        end
        
        
    case 'Ave'
        for i_src = 1 : size(src_loc, 2)
            dis(i_src) = mean(norms(repmat(src_loc(:, i_src), [1, size(dst_loc, 2)]) - dst_loc));            
        end
        
        
    case 'All'
        for i_src = 1 : size(src_loc, 2)
            dis(:, i_src) = norms(repmat(src_loc(:, i_src), [1, size(dst_loc, 2)]) - dst_loc);            
        end
        
        
    case 'Min+'
        for i_src = 1 : size(src_loc, 2)
            dummy  = norms(repmat(src_loc(:, i_src), [1, size(dst_loc, 2)]) - dst_loc);
            if(~isempty(dummy(dummy > 0)))
                dis(i_src)  = min(dummy(dummy > 0)); 
            else
                dis(i_src)  = 0;
            end
        end
       
        
    case 'Ave-Min+'
        dis    = 0;
        for i_src = 1 : size(src_loc, 2)
            dummy  = norms(repmat(src_loc(:, i_src), [1, size(dst_loc, 2)]) - dst_loc);
            if(~isempty(dummy(dummy > 0)))
                dis(i_src)  = min(dummy(dummy > 0)); 
            else
                dis(i_src)  = 0;
            end
        end
        dis  = mean(dis);
        
        
        
    case 'Ave-Min+-v2'
        count  = 0;
        dis    = 0;
        for i_src = 1 : size(src_loc, 2)
            dummy  = norms(repmat(src_loc(:, i_src), [1, size(dst_loc, 2)]) - dst_loc);
            count  = count + nnz(dummy > 0);
            dis    = dis + sum(dummy(dummy > 0));
        end
        dis  = dis / count;
        
    otherwise
        display('*************** Weird Request ***************')
        
        
end

end