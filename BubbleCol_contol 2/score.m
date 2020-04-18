function [total] = score(E,A)

    T1 = sum(E); % total EtOH produced
    
    T2 = 0; % penalize selectivity
 
    for i = 1:length(A)
        if A(i)<=.85 && A(i)>= .01 %penalize for .01<S<.85
            T2 = T2 + 1/A(i);
        elseif A(i) <= .01  %very stiff penalty for S<.01
            T2 = T2 + 1000; 
        end
    end
    
    total = T1-T2;
    
    
    
    
    
    