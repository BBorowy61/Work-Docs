function Vsh  = Vshift( V, S )
%   Vshift function shifts the vector elements,V, by the sepecified number
%   of positions, S down. 
%   The new vector last elements are padded with zeros
%Bogdan Borowy, EDF, 2016
Vsh = zeros(size(V));    
    for i=1:length(V),
        if ((i+S)>length(V) || (i+S)<1),
            Vsh(i) = 0;
        else Vsh(i)=V(i+S);
        end
    end
end

