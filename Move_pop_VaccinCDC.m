% dSi,a,r,vu =  -Poisson(Si,a,r,vui,a,r)
% dSi,a,r,vv =  Poisson(Si,a,r,vui,a,r)

% delayVac = 14;
delayT = t -delayVac*hourlyPerD;
if delayT>=1
    % first vaccin
    tmp_Pop = S(delayT, i, a, r, v)  + L(delayT, i, a, r, v)+ C(delayT, i, a, r, v)...
        + Ia(delayT, i, a, r, v)+Iy(delayT, i, a, r, v)+It(delayT, i, a, r, v) + R(delayT, i, a, r, v);

    % for S
    temp = V(delayT, i, a, r)*S(delayT, i,a,r,v)/tmp_Pop;
    
    if isnan(temp)
        temp = 0;
    end
    
    if S(delayT, i,a,r,v)>0 && temp>0
        S(t-1, i,a,r,v) = S(t-1, i,a,r,v) - temp*S(t-1, i,a,r,v)/ S(delayT, i,a,r,v);
        
        if   S(t-1, i,a,r,v)<0 || isnan(S(t-1, i,a,r,v)<0)
            S(t-1, i,a,r,v)=0;
        end
        
        S(t-1, i,a,r,v+1) = S(t-1, i,a,r,v+1) + temp*S(t-1, i,a,r,v)/ S(delayT, i,a,r,v);
    
    end
end
