function epi = revise_sdp(epi, sdi)
% revises seed position in the vector
% brings seed index to the head in entire patch indices vector 
for i = 1 : length(epi)
    epi{i}(epi{i} == sdi(i)) = [];
    epi{i} = [sdi(i) epi{i}];
end

end