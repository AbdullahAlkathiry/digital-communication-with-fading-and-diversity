% EE 419 - 212 - Term Project
% This project is the work of:
% Ahmed Alghamdi (Student A)
% Abdullah AlKathiry (Student B and C)
% Hussain Balhareth (Student C)

clear; clc; clf;
%% Simulation Settings ( only change this part of the code all needed result)
modelation = "MPSK"; % MPSK, QAM, or MPAM
M_modelation = 4;
recording = true; % Takes time to process
has_fading = true;
diversity_order = 4;% should be 1 if fading is not applied
div_type = 'MRC' ; % MRC or EGC
SNR_voice = 11;% if recording true, choose the SNR you want to hear voice at
%% producing the bits sequence
if recording
    fs=8000; %Sampling frequency
    recObj = audiorecorder(fs,24,1);
    disp('Start speaking.');
    recordblocking(recObj, 3);
    disp('End of Recording.');
    play(recObj);% Play back the recording
    Original = getaudiodata(recObj);% Store data in double-precision array.
    % convert to gray code binary quntized to
    signal = int16((2^9)*(1+(Original/max(abs(Original)))));
    x = reshape(dec2gc(signal,11)',1,[]);
    N_of_bits = length(x);
else
    N_of_bits = 1E4 + 128;
    x=round(rand(1,N_of_bits));
end
%% coordinates of the modulation 
if modelation == "MPSK"
    bits_per_symbol = log2(M_modelation);
    coordinate = zeros(M_modelation,2);
    phase = 2*pi/M_modelation;
    %Seting up the coordinates of MPSK
    for j = 0:M_modelation-1
        coordinate(j+1,:) = [cos(j*phase),sin(j*phase)];
    end
elseif modelation == "QAM"
    bits_per_symbol = log2(M_modelation);
    %Seting up the coordinates of QAM using the function qammod
    complex = qammod(0:M_modelation-1,M_modelation);
    coordinate = [real(complex'),imag(complex')];
    
elseif modelation == "MPAM"
    bits_per_symbol = log2(M_modelation);
    coordinate = zeros(M_modelation,2);
    %Seting up the coordinates of MPAM
    for j=1:M_modelation/2
        coordinate(j,:) = [2*j - 1 , 0];
        coordinate(j+ M_modelation/2 ,:) = [-2*j + 1 , 0];
    end
end
%% bits to symbols
w = rem(N_of_bits,bits_per_symbol); %Remainder
w1 = ceil(N_of_bits/bits_per_symbol);%number of symbols
if w~=0
    %We remove the extra bits if there are any
    x = [x,zeros(1,bits_per_symbol -w)];
end
N_of_bits = length(x);
y = zeros(w1,bits_per_symbol); %Setting up the symbols matrix
% bit gray coded to symbol
for j = 1:w1
    y(j,:)=x((j-1)*bits_per_symbol+1:j*bits_per_symbol);
end
s_num = gc2dec(y)+1;
coordinate_s =zeros(length(s_num),2);
for j = 1 : length(s_num)
    coordinate_s(j,:) = round(coordinate(s_num(j),:));
end
%% AWGN, Fading, and Diversity
Counter = 1;
for EbN0dB=1:12 %Vary SNR from 0 dB to 12dB
    EbN0=10^(EbN0dB/10); %SNR in linear
    N0=1/EbN0; %Calculate N0
    if has_fading
        %adding fading
        fading = sqrt(0.5)*randn([size(coordinate_s),diversity_order]);
        fading = sqrt(fading(:,1,:).^2 + fading(:,2,:).^2);
    else
        %not adding fading by making fading = 1
        fading = ones([size(coordinate_s),diversity_order]);
    end
    coordinate_f = fading.*coordinate_s; %Symbols on coordinate after fading
    %Adding AWGN
    noise = sqrt(N0/2)*randn(size(coordinate_f));
    coordinate_r = coordinate_f + noise; %Symbols on coordinate after fading + AWGN
    %Applying diversity
    if div_type == 'MRC'
        coordinate_r = cumsum(fading.*coordinate_r,3); %MRC
    elseif div_type == 'EGC'
        coordinate_r = cumsum(coordinate_r,3); %EGC
    end
    r_recived = zeros(length(s_num),diversity_order); %Matrix set  up
    %Calculating for multiple orders of diversity
    for i =1:diversity_order
        for j =1:length(s_num)
            [M,I]= min(sum((coordinate-coordinate_r(j,:,i)).^2,2));
            r_recived(j,i) = I ;
        end
    %Converting the symbols back to a bit stram
    bit_recived(i,:) = reshape(dec2gc(r_recived(:,i)-1,bits_per_symbol)',1,[]);
    end
    % Calculating SER and BER
    symbol_er(Counter,:) = sum((r_recived-s_num)~=0)/length(s_num); %SER
    bit_er(Counter,:) = sum((bit_recived-x)~=0,2)/N_of_bits; %BER
    EbN0dB_s(Counter)=EbN0dB; %The x-axis in the plot
    if recording && Counter == SNR_voice
       % To play the audio later on
       voice_recived = bit_recived(diversity_order,:);
    end
    clc;
    Counter=Counter+1
end
%% Plotting
clf
EbN0dB = -1:0.1:12; %logarithmic
EbN0=10.^(EbN0dB/10); %linear
if ~has_fading
    % Theoritical plots for AWGN channel
    d = sqrt(sum((coordinate(1,:)-coordinate(2,:)).^2));
    if modelation == "MPSK"
        Pe_= erfc(0.5*d*sqrt(EbN0));
        Pe_BER = Pe_/bits_per_symbol;
    elseif modelation == "QAM"
        Pe_= 2*(1-1/sqrt(M_modelation))*erfc(sqrt(EbN0));
    else
        Pe_= erfc(0.5*d*sqrt(EbN0));
        Pe_BER = Pe_/bits_per_symbol;
    end
    if M_modelation == 2
        semilogy(EbN0dB,Pe_BER,EbN0dB_s,bit_er,'x')
        legend('Theory BER','Simulation BER','Location','southwest');
    elseif modelation == "QAM"
        semilogy(EbN0dB,Pe_,EbN0dB_s,symbol_er,'o',EbN0dB_s,bit_er,'x')
        legend('Theory SER','Simulation SER','Simulation BER','Location','southwest');
    else
        semilogy(EbN0dB,Pe_,EbN0dB,Pe_BER,EbN0dB_s,symbol_er,'o',EbN0dB_s,bit_er,'x')
        legend('Theory SER','Theory BER','Simulation SER','Simulation BER','Location','southwest');
    end
    title(modelation + " with M = " +M_modelation + ", AWGN channel");
elseif M_modelation == 2 && modelation ~= "MPAM"
    % Theoritical plot of BER over fading channel for BPSK only
    Pe_fading = 0.5*(1-sqrt(EbN0./(1+EbN0)));
    semilogy(EbN0dB,Pe_fading); hold on; %Theory BER
    semilogy(EbN0dB_s,symbol_er,'x');%Sim BER
    legend(["Theory BER with div = 1" "Simulation BER with div ="+ [1:diversity_order]],'Location','southwest');
    title(modelation + " with M = " +M_modelation + ", fading channel, and "+ div_type +" diversity order from 1 to " + diversity_order);
else
    semilogy(EbN0dB_s,symbol_er,'o',EbN0dB_s,bit_er,'x')
    %legend("Fading SER", "Fading BER");
    legend(["SER with div ="+ [1:diversity_order] , "BER with div ="+ [1:diversity_order]],'Location','southwest');
    title(modelation + " with M = " +M_modelation + ", fading channel, and "+ div_type +" diversity order from 1 to " + diversity_order);
end
ylim([1E-6 1]); grid on;
xlabel("Eb/N0 (dB)"); ylabel("Probability of Error, Pe");
%% Play the recording 
if recording
    signal_recoverd = gc2dec(reshape(voice_recived,11,[])');
    signal_recoverd1 = (signal_recoverd/(2^9) - 1)*max(abs(Original));
    soundsc(signal_recoverd,fs)
end