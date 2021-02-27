function ang  = find_angle(src_ori, dst_ori, ang_opt)

ang  = zeros(size(src_ori, 2), size(dst_ori, 2));

for i = 1 : size(src_ori, 2)
    
    a  = src_ori(:, i);
    
    for j = 1 : size(dst_ori, 2)
        
        b  = dst_ori(:, j);
        ang(i, j)  = atan2(norm(cross(a,b)), dot(a,b))*180/pi;
        
    end
    
end

switch ang_opt
    case '[0, 180]'
        
    case '[0, 90]'
        ang(ang > 90)  = 180 - ang(ang > 90);
    
    case '[-90, 90]'
        ang(ang > 90)  = ang(ang > 90) - 180;
end

end



