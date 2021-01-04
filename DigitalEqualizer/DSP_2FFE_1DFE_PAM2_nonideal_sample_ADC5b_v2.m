figure(1);
clear all;clf;clc;
warning off;

% DSP Platform
%---------------------------------------
% Declare Parameter
    TMP_LEN     = 2^10-1;
    
    size_       = 10;

    Nu_of_FFE   = 2;
    Nu_of_DFE   = 1;
    
    muy         = 0.00000125;

% ADC no of bit and range of ADC   
    totalbits = 5;
    R         = [-31 31];
    
% Assume Non-Ideal sample frequency
%     scale = 4.5/10;
    scale     = 3/10;

% Assume parameter of channel
      g = [0.4 1 0.3]
      loss = 0.5;
      
%---------------------------------------

% PHY training sequence for PAM-2 Modulation
% ------------------------------------------------------------------------------
% Generate Test Mode 1 Per Gear, 2^21-2 sC2 balanced sequence
    sc4     = [15 -15];

% 2^20-1 symbol sequence used once at positive polarity and then in negative polarity to ensure zero DC
    scr     = ones(1,20);
    scrN    = zeros(TMP_LEN,20);

    for nn = 1:TMP_LEN
        scrN(nn,:) = scr;
        scr        = [mod(scr([1:16]+4)+scr([1:16]+1),2),scr(1:4)];
    end
    
    sc4Ai          = scrN(:,[3:4]+12);
    sc4Airef       = sc4Ai(:,1)+1;
    sc4Asymb       = sc4(sc4Airef);
    ideal_result   = sc4Airef'-1;
    inverse_ideal  = ones(1,size(ideal_result,2));
    
    for i = 1:size(ideal_result,2)
        if ideal_result(1,i) == 0
            inverse_ideal(1,i) = 1;
        else
            if ideal_result(1,i) == 1
                inverse_ideal(1,i) = 0;
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
    % Input signals
origin_input = input(1:size(input,2));

    % After Channel
after_channel = loss*filter(g,1,origin_input);

    % Non-ideal Sample
for i =1:size(input,2)-1
    nonidealsample(1,i) = (after_channel(1,i+1)-after_channel(1,i))*(scale+i)+after_channel(1,i)-(after_channel(1,i+1)-after_channel(1,i))*i;
end 

subplot(4,1,1);    
    stairs(origin_input,'r');
    oi = origin_input;
hold on;
grid on;
    plot(after_channel,'-*g');
    ac = after_channel;
    plot(nonidealsample,'-ob');
    ni = nonidealsample;
%------------------------------------------------------------------------------
% ADC 5 bit (Format: signed integer)
%------------------------------------------------------------------------------
X         = nonidealsample;
sz        = size(X);

if sz(1) > sz(2)
  X = X.';
end

q_lev = 2^totalbits -1;
midPnt = q_lev/2;    % center point
R_max = R(2);
R_min = R(1);
step = (R_max - R_min)/q_lev;
offset = 0.5*step;

% min_max clamping
X(find(X>=R(2))) = R_max;
X(find(X<=R(1))) = R_min;

% quantization
adc_out = (round((X-R_min)./step))*step;
adc_out = adc_out + R_min;
adc_out(find(adc_out > R_max-offset)) = R_max-step;
deltaR = diff(R);

if deltaR > 0  
  x_min = R(1) - deltaR/10;
  x_max = R(2) + deltaR/10;
else
  error('The upper limit must greater than the lower limit\n')
end

plot(adc_out,'-sk');
ao = adc_out;
axis([50 100 min(origin_input)-3 max(origin_input)+30]);

legend('Output of Tx','After Channel','After non-ideal sample','After ADC');

%-----------------------------------------------
    % Output of Tx
input = origin_input(1:size(origin_input,2)-100);
ideal_result = ideal_result(1:size(input,2));

    % Input of DSP
DSPin = adc_out(1:size(input,2));
tap = Nu_of_FFE + Nu_of_DFE;

    % Initial coefficient
wopt = zeros(tap,length(input)+1);
wopt(:,tap) = zeros(tap,1);
uvec = ones(tap,1);
delay1 = 0;

for i=tap:length(input)
    %   Bring new input sample
    ffe = DSPin(1,i)*wopt(2,i) + DSPin(1,i-1)*wopt(3,i);
    
    switch(delay1)
        case 0
            temp = ffe - 15*wopt(1,i);
        case 1
            temp = ffe + 15*wopt(1,i);
    end
        
    if (temp >= 0)
        test = 0;
    else 
        test = 1;
    end

    dfe = test;

    uvec(1) = delay1;
    uvec(2) = DSPin(1,i);
    uvec(3) = DSPin(1,i-1);

    %   Calculate the error
    e(i) = input(1,i-2) - temp;
    e2(i) = ideal_result(1,i-2)-dfe;
    
    delay1 = dfe;
    
    %   Update weights
    wopt(:,i+1)=wopt(:,i)+muy*uvec*e(i);

end


% Plot FFE and DFE parameter
    subplot(4,1,2);
    plot([1:length(input)+1],wopt(1,:),'b');
    hold on;
    grid on;
    plot([1:length(input)+1],wopt(2,:),'r');
    plot([1:length(input)+1],wopt(3,:),'g');
    down = min(min(wopt));
    up   = max(max(wopt));
    axis([1 length(input)+1 down-0.283 up+2.5]);
    legend('DFE Coefficient','FFE Coefficient (No Delay)','FFE Coefficient (1 Delay)');
    xlabel('Iteration');ylabel('Weights');
    title(['DSP Training Coefficient; Mu = ' num2str(muy)]);
    disp(['FFE coefficient: FFE = ', mat2str(wopt(2:3,length(input)+1))]);
    disp(['DFE coefficient: DFE = ', mat2str(wopt(1,length(input)+1))]);

    subplot(4,1,3);
    plot(e);
    hold on;
    grid on;
    plot(e2);
    legend('Error to modify LMS','Error between input sequence and DSP out');
    xlabel('Iteration');ylabel('Error');
    title('Error Estimation');
    axis([0 size(e,2) min(e)-1 max(e)+1]);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply result

    a           = wopt(1,length(input));
    b           = wopt(2,length(input));
    c           = wopt(3,length(input));
   
    N_START     = 16;
    N_END       = 20000;
    
    input       = input(N_START:N_END);
    ideal_result= ideal_result(N_START:N_END);

    after_channel = adc_out(N_START:N_END);
    
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
                temp = ffe(1,i) + 15*a;
        end

       if (temp >= 0)
            test = 0;
        else 
            test = 1;
        end

        dfe = test;
    
        delay1          = dfe;
        output          = [output dfe];
        before_decision = [before_decision temp];
    end
    
    subplot(4,1,4);
    stairs(ideal_result);
    hold on
    grid on
    stairs(output+2);
    axis([50 100 -0.5 5]);
    legend({'Transmit Sequence','Output + 2 Sequence'});
    xlabel('Sample Time');ylabel('Value');
    title(['Input and (Output + 2) Waveform']);

% Evaluate the error
No_error_bits    = sum(output(6:end)~=ideal_result(6:size(output,2)))
FER_in_percent   = No_error_bits/(size(output,2)-5)*100