function [E,rssum] = forcecurve(num, E, k, v,beta0)
R=8;
rssum = 0;
betae = [beta0,0];
parfor i=1:num*num
    FileName1 = strcat('Line',num2str(floor((i-1)/num),'%04d'),'Point',num2str(i-(floor((i-1)/num))*num-1,'%04d'),'.ibw');
    a=IBWread(FileName1);
    raw=a.y(:,1);
    defl=a.y(:,2)*10^9;
    zsnsr=a.y(:,3)*10^9;
    
    %         force=a.y(:,4)*10^9;
    %         ind=a.y(:,5)*10^9;
    
    force = defl*k;
    ind = zsnsr-defl;

    % plot(ind,force);
    % xlabel('Indentation/nm');
    % ylabel('Force/nN');

    
    [~,maxindex]=max(zsnsr);
    
    afor = force(1:maxindex);    % approach
    aind = ind(1:maxindex);
    rfor = force(maxindex+1:end);% retract
    rind = ind(maxindex+1:end);
    rfor = rfor(end:-1:1);
    rind = rind(end:-1:1);
    
    % plot(aind,afor);
    % hold on
    % plot(rind,rfor);
    % xlabel('Indentation/nm');
    % ylabel('Force/nN');
    % legend('approach curve','retract curve','Location','northwest');
    % hold off

    %         [m,n]=min(a(:,2));
    %         MIN=a(n,:);
    zero_force_range = [aind(round(0.2*length(aind)):round(0.6*length(aind))),afor(round(0.2*length(aind)):round(0.6*length(aind)))];
    zero_force_fit = polyfit(zero_force_range(:,1),zero_force_range(:,2),1);
    azero_force_eval = zero_force_fit(1).*aind + zero_force_fit(2);
    rzero_force_eval = zero_force_fit(1).*rind + zero_force_fit(2);
    temzafor = afor - azero_force_eval;
    temzrfor = rfor - rzero_force_eval;

    % plot(aind,temzafor);
    % hold on
    % plot(rind,temzrfor);
    % xlabel('Indentation/nm');
    % ylabel('Force/nN');
    % legend('preliminarily zeroized approach curve','preliminarily zeroized retract curve','Location','northwest');
    % hold off

    stdev1 = std(afor(1:round(0.3*length(aind))));
    
    %fType 1 means retract falls first, 2 means approach falls first
    [asp,rsp,fType,gof] = sinfit(temzafor,temzrfor,aind,rind,stdev1);
    if gof.rsquare < 0.4
        gof.rsquare;
    end
    rssum = rssum+gof.rsquare;

    switch fType
        case 1
            [~,asp] = min(abs(aind-rind(rsp)));
        case 2
            [~,rsp] = min(abs(rind-aind(asp)));
    end
    zrfor = [temzrfor(1:rsp);rfor(rsp+1:end)-rfor(rsp)];
    zafor = [temzafor(1:asp);afor(asp+1:end)-afor(asp)];
    [minafor,minaforindex] = min(zafor(asp+1:end));
    [~,minrforindex] = min(zrfor);
   
    if minafor<-2*gof.rmse
        for m = minaforindex+asp-1:length(zafor)
            if zafor(m)<0 && zafor(m+1)>=0
                CP=m;
                break;
            end
        end
        
        % if no ditch fit to straight line
    else
        CP=9999;
        fullLen = length(aind)-asp;
        fitRangeSt = round(asp+0.05*fullLen); % start
        fitRangeEn = round(asp+0.3*fullLen); % end
        stfit = polyfit(aind(fitRangeSt:fitRangeEn),zafor(fitRangeSt:fitRangeEn),1);
        k2 = stfit(1);
        b2 = stfit(2);
        for m = asp-10:length(zafor)
            if aind(m-1) + b2/k2 < 0 && aind(m) + b2/k2 >=0
                CP = m;
                break;
            end
        end
        if CP == 9999
            CP = asp;
        end
    end

    zeroind = (aind(CP)+aind(CP+1))/2;
    zaind = aind - zeroind;
    YFitLen = length(zaind(CP:end));
    YFitSt = round(CP+0.1*YFitLen);
    YFitEn = round(CP+0.5*YFitLen);
    YFitInd = zaind(YFitSt:YFitEn);
    YFitFor = zafor(YFitSt:YFitEn);

    Fad = abs(min(zrfor));
    %% Hertz model
    % fun = @(x,fforce) R^(2/3)*(fforce+Fad).^(2/3)/x(1)^(2/3)/R+x(2);

    %% DMT model
    fun = @(x,fforce) R^(2/3)*fforce.^(2/3)/x(1)^(2/3)/R+x(2);

    etot=nlinfit(YFitFor,YFitInd,fun,betae);
    e = (1-v^2)/(4/3/etot(1)-(1-0.25^2)/(169));
    
    E(i)=e;
    fclose all;
end