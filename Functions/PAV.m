function Pmono=PAV(PP)

    n1 = max(size(PP));
    wt = ones(size(PP));
    lvlsets = (1:n1)';
    one2 = (1:n1-1)';

    while (true)
        viol = (PP(1:(n1-1))-PP(2:n1) > 0);
        if (~any(viol))
            break;
        end
        ii = min(find(viol));
        lvl1 = lvlsets(ii);
        lvl2 = lvlsets(ii+1);
        ilvl = (lvlsets==lvl1)|(lvlsets==lvl2);
        PP(ilvl) = sum( PP(ilvl).*wt(ilvl) )/sum(wt(ilvl));
        lvlsets(ilvl) = lvl1;
    end
    
    Pmono=PP;
    
    
    
 

   
