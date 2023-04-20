clc;
close all;
warning('off');
% E=zeros(1,10);
   %（Tip radius）
k=61.73; %spring constant Al=28.234 location1=64.554 location2=61.73 210309=3.114
           % poisson ratio from others
vs = 0.38;
beta0 = 1;
error = 0;
filepath = 'Epoxy/location_2';
addpath(filepath);
addpath('Igor2Matlab');
file = dir(filepath);
filenumber = size(file)-2;
num = sqrt(filenumber(1));
E = zeros(num*num,1);

[E,rssum] = forcecurve(num, E, k,vs,beta0);

for i=1:num*num
    if ~(isreal(E(i)))
        E(i) = E(i-1);
        error = error+1;
    elseif E(i) < 0
        E(i) = E(i-1);
        error = error+1;
    end
end
if error > 0.1*num*num
    disp('Too many errors in this dataset, debug it');
    exit;
end

E=E*10^9;
me = mean(E);
if me>10^9
    E = E/10^9;
    unit = 'GPa';
elseif me>10^6 && me<10^9
    E = E/10^6;
    unit = 'MPa';
elseif me>10^3 && me<10^6
    E = E/10^3;
    unit = 'kPa';
else
    E = E;
    unit = 'Pa';
end
box = E;
youngs = figure(1);
title('Young''s modulus');
plot(E);
xlabel('number');
lb = ['Young''s modulus/',unit];
ylabel(lb);
E=reshape(E,num,num);
E=rot90(E);
heatmap = figure(2);
imagesc(E);
pic=colorbar;
set(get(pic,'Title'),'string',unit);
colormap(gray);
title('Young''s Modulus');
xlabel('point');
ylabel('point');
bp = figure(3);
boxplot(box);

saveas(youngs,'Youngsmodulus_el2','fig');
saveas(heatmap,'heatmap_el2','fig');
saveas(bp,'boxplot_of_el2','fig');
save('boxel2.mat','box')
exit;
