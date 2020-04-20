function [total, viol] = score(E,A)

    total = sum(E); % total EtOH produced
    
    viol = 0; % penalize selectivity
 
    for i = 1:length(A)
        if E(i)/A(i)<=.85 
            viol = viol + 1;
        end
    end
    
end
    
    
    
    
    
    