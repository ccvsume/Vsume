function sp = sepoint(ind,smtemzfor,gof,pr)
fitfor = pr(ind);
for m = round(0.1*length(smtemzfor)):length(smtemzfor)
    
    if smtemzfor(m)<fitfor(m)-5*gof.rmse
        sp = m-1; % retract sp index
        break;
    elseif smtemzfor(m)>fitfor(m)+5*gof.rmse
        sp = m-1;
        break;
    end
end