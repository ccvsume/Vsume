function [asp,rsp,fType,gof] = sinfit(temzafor,temzrfor,aind,rind,stdev)

    [~,MinRForIndex] = min(temzrfor);
    [~,MinAForIndex] = min(abs(aind-rind(MinRForIndex)));
    if temzafor(MinAForIndex) < 3.5*stdev
        [rsmfit,~] = sm(rind,temzrfor);
        smtemzrfor = rsmfit(rind);
        [rpr,gof] = rsf(rind,temzrfor);
        rsp = sepoint(rind,smtemzrfor,gof,rpr);
        fType = 1;
        asp = 0;
    else
        [asmfit,~] = sm(aind,temzafor);
        smtemzafor = asmfit(aind);
        [apr,gof] = asf(aind,smtemzafor);
        asp = sepoint(aind,smtemzafor,gof,apr);
        fType = 2;
        rsp = 0;
    end

    % switch fType
    %     case 1
    %         plot(rind,temzrfor);
    %         hold on
    %         plot(rind,smtemzrfor);
    %         xlabel('Indentation/nm');
    %         ylabel('Force/nN');
    %         legend('preliminarily zeroized retract curve','smoothed preliminarily zeroized retract curve','Location','northwest');
    %         hold off
    %     case 2
    %         plot(aind,temzafor);
    %         hold on
    %         plot(aind,smtemzafor);
    %         xlabel('Indentation/nm');
    %         ylabel('Force/nN');
    %         legend('preliminarily zeroized approach curve','smoothed preliminarily zeroized approach curve','Location','northwest');
    %         hold off
    % end
end

function [pr, gof] = asf(indentation,force)

    [~,range] = min(force);
    if range<0.8*length(force)
        [~,range] = min(force(round(0.8*length(force)):end));
        range = range + floor(0.8*length(force));
        fitrange = findnearest(force(range-20:range));
        if fitrange == 0
            fitrange = range;
        else
            fitrange = range-20+fitrange;
        end
    else
        fitrange = findnearest(force(range-20:range));
        if fitrange == 0
            fitrange = range;
        else
            fitrange = range-20+fitrange;
        end
    end
    
    fitx = indentation(1:fitrange);
    fity = force(1:fitrange);
    [xData, yData] = prepareCurveData( fitx, fity );
    
    % 设置 fittype 和选项。
    ft = fittype( 'sin8' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf];
    % opts.StartPoint = [0.0198014770467271 0.0153116423949833 -2.10139824027484 0.0113737215500289 0.00765582119749165 2.18385723409194 0.00700777322988467 0.0114837317962375 1.32002170419094 0.00353429967019636 0.00191395529937291 2.86611482632828 0.0049513369249967 0.0306232847899666 2.32846763427691];
    % opts.StartPoint = [0.327121682384547 0.0139108926175445 3.10056375942429 0.243130768540258 0.00695544630877223 -1.86178470757003 0.0149067603338926 0.0208663389263167 2.47139675297498 0.862745259299006 0.00347772315438611 -0.912048663940849 0.00621969078026599 0.0278217852350889 2.93781544529696];

    % 对数据进行模型拟合。
    % 5point smooth数据拟合效果不如不smooth的
    [pr, gof] = fit( xData, yData, ft, opts );
end


function [pr, gof] = rsf(indentation,force)
    [~,range] = min(force);
    fitrange = findnearest(force(range-150:range))+range-150;
    
    fitx = indentation(1:fitrange);
    fity = force(1:fitrange);
    [xData, yData] = prepareCurveData( fitx, fity );
    
    % 设置 fittype 和选项。
    ft = fittype( 'sin5' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf];
    % opts.StartPoint = [0.0198014770467271 0.0153116423949833 -2.10139824027484 0.0113737215500289 0.00765582119749165 2.18385723409194 0.00700777322988467 0.0114837317962375 1.32002170419094 0.00353429967019636 0.00191395529937291 2.86611482632828 0.0049513369249967 0.0306232847899666 2.32846763427691];
    % opts.StartPoint = [0.327121682384547 0.0139108926175445 3.10056375942429 0.243130768540258 0.00695544630877223 -1.86178470757003 0.0149067603338926 0.0208663389263167 2.47139675297498 0.862745259299006 0.00347772315438611 -0.912048663940849 0.00621969078026599 0.0278217852350889 2.93781544529696];
   
    % 对数据进行模型拟合。
    % smooth数据拟合效果不如不smooth的
    [pr, gof] = fit( xData, yData, ft, opts );
end

function fitrange = findnearest(force)
    if isempty(force)
        fitrange = 0;
    elseif force(end-1) - force(end) > 0 || force(end-1) < -20
        force = force(1:end-1);
        fitrange = findnearest(force);
    else
        fitrange = length(force);
    end
    
end
