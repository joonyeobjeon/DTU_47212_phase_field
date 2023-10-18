function [matrix] = make_multi_phase(AA)
    NP = unique(AA);
    matrix = {};
    for i=1:size(NP,1)
        matrix=[matrix; {zeros(size(AA,1), size(AA,2))}];
    end

    for i=1:size(NP,1)
        A = find(AA==i)
        matrix{1}(:,:)=0
        matrix{1}(A)=1
    end

return