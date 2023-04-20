clc;
close all;
warning('off');
% E=zeros(1,10);
   %（Tip radius）
k=28.234; %spring constant Al=28.234 location1=64.554 location2=61.73 210309=3.114
           % poisson ratio from others
count = 0;
filepath = 'Al/Aluminium3';
addpath(filepath);
addpath('Igor2Matlab');
file = dir(filepath);
filenumber = size(file)-2;
num = sqrt(filenumber(1));
E = zeros(num*num,1);

[E,rssum] = forcecurve(num, E, k);

for i=1:num*num
    if ~(isreal(E(i)))
        E(i) = E(i-1);
        count=count+1;
    elseif i == 44*num+20
        E(i) = E(i-1);
        continue;
    end
end

E=E*10^9;
box = E;
youngs = figure(1);
title('E')
plot(E);
xlabel('number');
ylabel('E');
E=reshape(E,num,num);
E=rot90(E);
heatmap = figure(2);
imagesc(E);
pic=colorbar;
set(get(pic,'Title'),'string','Pa');
colormap(gray);
title('Young''s Modulus');
xlabel('point');
ylabel('point');
bp = figure(3);
boxplot(box);

saveas(youngs,'Youngsmodulus_Al','fig');
saveas(heatmap,'heatmap_Al','fig');
saveas(bp,'boxplot_of_Al','fig');

% exit;
