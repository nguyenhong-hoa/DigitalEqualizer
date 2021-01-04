%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Digital Signal Processing finding FFE & DFE co-eff
% Scan parameter is placed in v1
% Sample Pattern (Suggestted by A Phi) is placed in v1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;clf;clc;
warning off;
% DSP Platform
%---------------------------------------
    TMP_LEN     = 2^10-1;
    
    size_       = 12;
    
    Nu_of_FFE   = 2;
    Nu_of_DFE   = 1;
    
    muy         = 0.00000125;

% Assume parameter of channel

%     g = [0.1 1 0.7]
%     g = [0.2 1 0.8]
%     g = [0.3 1 0.1]
%     g = [0.3 1 0.2]
%     g = [0.3 1 0.7]
%     g = [0.4 1 0.3]
    g = [0.4 1 0.5]

% PHY training sequence for PAM-4 Modulation
% ------------------------------------------------------------------------------
% Generate Test Mode 1 Per Gear, 2^21-2 sC2 balanced sequence
    sc4     = [15 5 -15 -5];

% 2^20-1 symbol sequence used once at positive polarity and then in negative polarity to ensure zero DC
    scr     = ones(1,20);
    scrN    = zeros(TMP_LEN,20);

    for nn = 1:TMP_LEN
        scrN(nn,:) = scr;
        scr        = [mod(scr([1:16]+4)+scr([1:16]+1),2),scr(1:4)];
    end
    
    sc4Ai          = scrN(:,[3:4]+12);
    sc4Airef       = (sc4Ai(:,1)*2+sc4Ai(:,2))+1;
    sc4Asymb       = sc4(sc4Airef);
    ideal_result   = sc4Airef'-1;
    inverse_ideal  = ones(1,size(ideal_result,2));
    
    for i = 1:size(ideal_result,2)
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

    Aseq         = [sc4Asymb -sc4Asymb];
    input        = Aseq;
    
    for i = 1:size_
        input = [input input];
        ideal_result = [ideal_result ideal_result];
    end
    
%------------------------------------------------------------------------------

% Extract Input Sequence
    input = input(1:size(ideal_result,2));
    ideal_result = ideal_result(1:size(input,2));
    after_channel = filter(g,1,input);

    tap = Nu_of_FFE + Nu_of_DFE;

% Initial coefficient
wopt = zeros(tap,length(input)+1);
wopt(:,tap) = zeros(tap,1);
uvec = ones(tap,1);
delay1 = 0;

for i=tap:length(input)
    %   Bring new input sample
    ffe = after_channel(1,i)*wopt(2,i) + after_channel(1,i-1)*wopt(3,i);
    
    switch(delay1)
        case 0
            temp = ffe - 15*wopt(1,i);
        case 1
            temp = ffe - 5*wopt(1,i);
        case 2
            temp = ffe + 15*wopt(1,i);
        case 3
            temp = ffe + 5*wopt(1,i);
    end
        
    v1 = 0;
    v2 = 0;
    v3 = 0;

    if (temp >= 10)
        v1 = 1;
    end

    if (temp >= 0)
        v2 = 1;
    end

    if (temp >= -10)
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

    uvec(1) = delay1;
    uvec(2) = after_channel(1,i);
    uvec(3) = after_channel(1,i-1);

    %   Calculate the error
    e(i) = input(1,i-2) - temp;
    e2(i) = ideal_result(1,i-2)-dfe;
    
    delay1 = dfe;
    
    %   Update weights
    wopt(:,i+1)=wopt(:,i)+muy*uvec*e(i);

end

% Plot FFE and DFE parameter
    figure(1);
    clf;
    subplot(2,1,1);
    plot([1:length(input)+1],wopt(1,:),'b');
    hold on;
    grid on;
    plot([1:length(input)+1],wopt(2,:),'r');
    plot([1:length(input)+1],wopt(3,:),'g');
    down = min(min(wopt));
    up   = max(max(wopt));
    axis([1 length(input)+1 down-0.283 up+0.877]);
    legend('DFE Coefficient','FFE Coefficient (No Delay)','FFE Coefficient (1 Delay)');
    xlabel('Iteration');ylabel('Weights');
    title(['DSP Training Coefficient; Mu = ' num2str(muy)]);
    disp(['FFE coefficient: FFE = ', mat2str(wopt(2:3,length(input)+1))]);
    disp(['DFE coefficient: DFE = ', mat2str(wopt(1,length(input)+1))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
% Apply result

    a           = wopt(1,length(input));
    b           = wopt(2,length(input));
    c           = wopt(3,length(input));
   
    N_START     = 20;
    N_END       = 20000;
    
    input       = input(N_START:N_END);
    ideal_result= ideal_result(N_START:N_END);

    after_channel = filter(g,1,input);
    
    for i=4:size(input,2)
        ffe(1,i)  = after_channel(1,i)*b + after_channel(1,i-1)*c;
    end
    
    output      = [];
    before_decision = [];
    delay1      = 0;
    
    for i=3:size(input,2)
        switch(delay1)
            case 0
                temp = ffe(1,i) - 15*a;
            case 1
                temp = ffe(1,i) - 5*a;
            case 2
                temp = ffe(1,i) + 15*a;
            case 3
                temp = ffe(1,i) + 5*a;
        end
        
        v1  = 0;
        v2  = 0;
        v3  = 0;
        
        if (temp >= 10)
            v1 = 1;
        end
        
        if (temp >= 0)
            v2 = 1;
        end
        
        if (temp >= -10)
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
        delay1          = dfe;
        output          = [output dfe];
        before_decision = [before_decision temp];
    end
    
subplot(2,1,2);
stairs(ideal_result);
hold on
grid on
stairs(output+4);
axis([50 100 -0.5 9.5]);
legend({'Transmit Sequence','Output + 4 Sequence'});
xlabel('Sample Time');ylabel('Value');
title(['Input and (Output + 4) Waveform']);

% Evaluate the error
No_error_bits = sum(output(6:end)~=ideal_result(6:size(output,2)))
FER_in_percent   = No_error_bits/(size(output,2)-5)*100



