clear all
figure
hold on
colors = ['r','b','m','c','y'];
for i=0:4
    fileName = sprintf('output_%d.txt',i);
    content = fileread(fileName);
    data = textscan(content, '%f %f %f', 'HeaderLines', 0);
    
    t = data{1};
    v1 = data{2};
    v2 = data{3};
    
    % for validation
    %%time1 = linspace(0,7,2000000);
    %%realX = zeros(1,length(time1));
    %have t want to generate plot of actual solution of x
    %%for k=1:length(time1)
        %%time = time1(k);
        %%realX(k) = (4/1.3)*(exp(0.8*time)-exp(-0.5*time))+2*exp(-0.5*time);
    %%end
    
    plot(t,v1,colors(i+1));
    plot(t,v2, colors(i+1));
end
%legend('V1','V2');
xlabel('Time');
ylabel('Output');
title('Output vs. Time for Validation');
hold off
