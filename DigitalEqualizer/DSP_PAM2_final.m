%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Digital Signal Processing finding FFE & DFE co-eff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;clf;clc;
warning off;
% DSP Platform
%---------------------------------------
% Assume parameter of channel
g = [0.5 1 0.5]

% Pseudo Sequence Input
% Number of bit:
N = 5000;
% Option 1: sequence test clock
% test = [0 0 1 1];
% input =[];
% for i = 1:round(N/4)
%     input = [input test];
% end

% Option 2: sequence test random
input = randi([0 1],1,N);

% Option 3: sequence test PRBS
% O = 4;
% N = 2^round(log2(N))-1;
% input = prbs(O,N);

% Pseudo Output of Channel
output = filter(g,1,input);

% Ideal coefficient
a = -g(1,1);
b = g(1,2);
c = g(1,3);

% Initial FFE result
ffe = [a*output(1,1),zeros(1,(size(output,2)-1))];

% FFE
for i = 2 :size(output,2)
    ffe(1,i) = a*output(1,i) + b*output(1,i-1);
end

% Intial DFE result
dfe = zeros(1,size(ffe,2));

% DFE
for i =2 :size(output,2)
    if (dfe(1,i-1) == 1)
        dfe(1,i) = ((ffe(1,i) - c)>0.1);
    else
        dfe(1,i) = ((ffe(1,i)+c)>(0+c));
    end
end

output = dfe;

% Evaluate result
difference_ideal = sum(input(1,1:(size(input,2)-2))~=output(1,3:size(input,2)))

%---------------------------------------

% Extract Input Sequence
input = input(1,1:(size(input,2)));
after_channel = filter(g,1,input);

% Initial coefficient
wopt = zeros(size(g,2),length(input)+1);
wopt(:,3)=[0.1;0.1;0.1];
uvec = zeros(size(g,2),1);
muy = 1/2^8;

previous_bit = 0;

for i=3:length(input)
    % Bring new input sample
    ffe = wopt(1,i)*after_channel(1,i) + wopt(2,i)*after_channel(1,i-1);

        if (previous_bit == 1)
            dfe = ((ffe - wopt(3,i-1))>0.1);
        else
            dfe = ((ffe + wopt(3,i-1))>(0+wopt(3,i-1)));
        end
        
    previous_bit = dfe;
    
    uvec(1) = dfe;
    
    %Calculate the error
    e(i) = input(i-2)-uvec(1);
    
    %Update weights
    wopt(:,i+1)=wopt(:,i)+muy*uvec*conj(e(i));
    
    %shift the data sample
    uvec(3) = uvec(2);
    uvec(2) = uvec(1);
end

% Plot w parameter optimal vector versus Iteration
figure(1);
hold off;
plot([1:length(input)+1],wopt(1,:),'b');
hold on;
plot([1:length(input)+1],wopt(2,:),'r');
plot([1:length(input)+1],wopt(3,:),'g');
legend('Weight vector (1)','Weight vector (2)','Weight vector (3)');
xlabel('Iteration');ylabel('Weights');
grid on;
title(['LMS algorithm; Mu = ' num2str(muy)]);
disp(['Weight vector of AF: w = ', mat2str(wopt(:,length(input)+1))]);




% Apply result
a = wopt(1,length(input)+1);
b = wopt(2,length(input)+1);
c = wopt(3,length(input)+1);

input = randi([0 1],1,2000000);
output = filter(g,1,input);

ffe = [a*output(1,1),zeros(1,(size(output,2)-1))];

for i = 2 :size(output,2)
    ffe(1,i) = a*output(1,i) + b*output(1,i-1);
end

dfe = zeros(1,size(ffe,2));

for i =2 :size(output,2)
    if (dfe(1,i-1) == 1)
        dfe(1,i) = ((ffe(1,i) - c)>0.1);
    else
        dfe(1,i) = ((ffe(1,i) + c)>((0)+c));
    end
end

output = dfe;

difference_new = sum(input(1,1:(size(input,2)-2))~=output(1,3:size(input,2)))
