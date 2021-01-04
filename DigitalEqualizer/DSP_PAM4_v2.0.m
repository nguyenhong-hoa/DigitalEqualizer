%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Digital Signal Processing finding FFE & DFE co-eff
% Scan parameter is placed in v1
% Sample Pattern (Suggestted by A Phi) is placed in v1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;clf;clc;
warning off;
% DSP Platform
%---------------------------------------
    TMP_LEN=2^10-1;
    
    size_ = 5; %size = 20*(2^10-1)
    
    Nu_of_FFE = 3;
    Nu_of_DFE = 1;
    
    thres1 = 15;
    thres2 = 6;
    thres3 = 9;
    
    muy =0.000009;

% Assume parameter of channel

    % GOOD: 0%
    % g = [0.1 1 0.7]
    % g = [0.2 1 0.5]
    % g = [0.3 1 0.7]
    % g = [0.3 1 0.2]

    % 2.9%
    g = [0.3 1 0.1]

    % 3.5%
    % g = [0.4 1 0.5]

    % 6.3%
    % g = [0.4 1 0.3]

    % 11%
    % g = [0.4 1 0.2]

% PHY training sequence for PAM-4 Modulation
% ------------------------------------------------------------------------------
% Generate Test Mode 1 Per Gear, 2^21-2 sC2 balanced sequence
    sc4=[15 5 -15 -5];

% 2^20-1 symbol sequence used once at positive polarity and then in negative polarity to ensure zero DC
    scr=ones(1,20);
    scrN=zeros(TMP_LEN,20);

    for nn=1:TMP_LEN
        scrN(nn,:)=scr;
        scr = [mod(scr([1:16]+4)+scr([1:16]+1),2),scr(1:4)];
    end
    
    sc4Ai = scrN(:,[3:4]+12);
    sc4Airef = (sc4Ai(:,1)*2+sc4Ai(:,2))+1;
    sc4Asymb = sc4(sc4Airef);
    ideal_result = sc4Airef'-1;
    inverse_ideal = ones(1,size(ideal_result,2)); %No DC
    
    for i=1:size(ideal_result,2)
        if ideal_result(1,i) == 0
            inverse_ideal(1,i) = 2;
        else
            if ideal_result(1,i) == 1
                inverse_ideal(1,i) = 3;
            else
                if ideal_result(1,i)==2
                    inverse_ideal(1,i)=0;
                else
                    inverse_ideal(1,i)=1;
                end
            end
        end
    end
    
    ideal_result = [ideal_result inverse_ideal];

    Aseq = [sc4Asymb -sc4Asymb];
    input = Aseq;
    
    for i =1:size_
        input = [input input];
        ideal_result =[ideal_result ideal_result];
    end
    
%------------------------------------------------------------------------------

% Extract Input Sequence
    input = input(1:size(ideal_result,2)-200);

    after_channel = filter(g,1,input);

    tap = Nu_of_FFE + Nu_of_DFE;

% Initial coefficient
wopt = zeros(tap,length(input)+1);
wopt(:,3)=0.5*ones(tap,1);
uvec = ones(tap,1);

for i=tap:length(input)
    % Bring new input sample
    ffe = after_channel(1,i)*wopt(1,i) + after_channel(1,i-1)*wopt(2,i) + after_channel(1,i-2)*wopt(3,i) + after_channel(1,i-3)*wopt(4,i);
    
    uvec(1) = ffe;

    %Calculate the error
    e(i) = (input(1,i-2)-uvec(1));
    
    %Update weights
    wopt(:,i+1)=wopt(:,i)+muy*uvec*e(i);
    
    %shift the data sample
    uvec(3) = uvec(2);
    uvec(2) = uvec(1);
end

% Plot w parameter optimal vector versus Iteration
%     figure(1);
    subplot(2,1,1);
    hold off;
    plot([1:length(input)+1],wopt(1,:),'b');
    hold on;
    grid on;
    plot([1:length(input)+1],wopt(2,:),'r');
    plot([1:length(input)+1],wopt(3,:),'g');
    plot([1:length(input)+1],wopt(4,:),'k');
    legend('Weight vector (1)','Weight vector (2)','Weight vector (3)');
    xlabel('Iteration');ylabel('Weights');
    title(['LMS algorithm; Mu = ' num2str(muy)]);
    disp(['Weight vector of AF: w = ', mat2str(wopt(:,length(input)+1))]);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
% Apply result

    a = floor(wopt(2,length(input))*2^8)/2^8;
    b = floor(wopt(3,length(input))*2^8)/2^8;
    c = floor(wopt(1,length(input))*2^8)/2^8;
    d = floor(wopt(4,length(input))*2^8)/2^8;
%     c = abs(wopt(1,length(input))) + abs(wopt(2,length(input))) + abs(wopt(3,length(input))) 
%     a = wopt(3,length(input))
%     b = -wopt(4,length(input))
%     d = sum(wopt(:,length(input)))
    
    N_START = 1000;
    N_END   = 5000;
    
    input=input(N_START:N_END);
    ideal_result=ideal_result(N_START:N_END);

    after_channel = filter(g,1,input);
    
    for i=3:size(input,2)
        ffe(1,i) = after_channel(1,i-1)*a+after_channel(1,i)*b + after_channel(1,i-2)*c;
    end
    
    output = [];
    before_decision = [];
    delay1=0;
    
    for i=2:size(input,2)
        switch(delay1)
            case 0
                temp = ffe(1,i) - thres1*d;
            case 1
                temp = ffe(1,i) - thres2*d;
            case 2
                temp = ffe(1,i) + thres1*d;
            case 3
                temp = ffe(1,i) + thres2*d;
        end
        
        v1 =0;
        v2 =0;
        v3 =0;
        
        if (temp >= thres3)
            v1 = 1;
        end
        
        if (temp >= 0)
            v2 = 1;
        end
        
        if (temp >= -thres3)
            v3 = 1;
        end
        
        test = 4*v1+2*v2+v3;
        switch (test)
            case 7
                dfe = 0;
            case 3
                dfe = 1;
            case 1
                dfe = 3;
            case 0
                dfe = 2;
        end
        delay1 = dfe;
        output = [output dfe];
        before_decision = [before_decision temp];
    end
    
% figure(2)
% clf;
subplot(2,1,2);
hold on
grid on
% subplot(3,1,1);
stairs(ideal_result);
stairs(output+4);
axis([1 50 -1 8]);
legend('Transmit Sequence','Output Sequence');
xlabel('Sample Time');ylabel('Value');
title(['Input and (Output + 4) Waveform']);
% subplot(3,1,2);
% stairs(output);
% axis([0 30 -2 7]);
% grid on

% subplot(3,1,3);
% stairs(ideal_result);
% hold on
% grid on
% axis([0 20 -15 15]);
% stairs(before_decision);
 
error = sum(output(5:end)~=ideal_result(4:end-2))
FER   = error/size(input,2)

