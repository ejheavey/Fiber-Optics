%% Analyzing BER vs SNR by adjusting channel gain

M256 = 256;
k256 = log2(M256);
n = 30000;
samplesPerSymbol256 = 1;

rng default;

input_data256 = randi([0 1], n, 1);
input_data_matrix256 = reshape(input_data256, length(input_data256)/k256, k256); % reshape the data to n x n 
input_data_symbols256 = bi2de(input_data_matrix256); % binary to decimal conversion
% figure; stem(input_data_symbols256(1:10)); % show the input symbols in decimal

dataMod = qammod(input_data_symbols, M256, 'bin'); % modulate the input data using binary encoding
dataMod_GrayCode256 = qammod(input_data_symbols256, M256, 'UnitAveragePower', false); % ... and Gray coded
avg_mod_pwr256 = sum(abs(dataMod_GrayCode256).^2)/M256;
avg_mod_pwr256 = sum(abs(dataMod_GrayCode256).^2)/(length(dataMod_GrayCode256)/M256);

EbN0256 = 10;
snr0_256 = EbN0256 + 10*log10(k256) - 10*log10(samplesPerSymbol256);
receivedSignal_GC256 = awgn(dataMod_GrayCode256, snr0_256, 'measured');
s_plot = scatterplot(receivedSignal_GC256, 1, 0, 'g.');
hold on
scatterplot(dataMod_GrayCode256, 1, 0, 'k*', s_plot); title('Constellation Diagram with AWGN for 256-QAM Signal')



EbN0256 = 0;
for i = 1:50
    snr256(i) = EbN0256 + 10*log10(k256) - 10*log10(samplesPerSymbol256);
    receivedSignal_GrayCode256 = awgn(dataMod_GrayCode256, snr256(i), 'measured');
    avg_receivedPwr256(i) = sum(abs(receivedSignal_GrayCode256).^2)/M256;
    output_data_symbols_GrayCode256 = qamdemod(receivedSignal_GrayCode256, M256);

    output_data_matrix_GrayCode256 = de2bi(output_data_symbols_GrayCode256, k256);
    output_data_GrayCode256 = output_data_matrix_GrayCode256(:);

    [numbErrors_GrayCode256, BER_GrayCode256] = biterr(input_data256, output_data_GrayCode256);
    logber256(i) = log2(BER_GrayCode256);
    ber256(i) = BER_GrayCode256;
    EbN0256 = EbN0256 + 1;
end


figure, plot(snr256(1:21), logber256(1:21)); xlabel('SNR (dB)'), ylabel('log2(BER) (dB)'), title('log(BER) vs SNR, 16-, 64- and 256-QAM')

%% Analyzing BER vs SNR by adjusting samples

M256 = 256;
k256 = log2(M256);
n = 30000;
samplesPerSymbol256 = 1;

rng default;

input_data256 = randi([0 1], n, 1);
input_data_matrix256 = reshape(input_data256, length(input_data256)/k256, k256); % reshape the data to n x n 
input_data_symbols256 = bi2de(input_data_matrix256); % binary to decimal conversion
figure; stem(input_data_symbols256(1:10)); % show the input symbols in decimal

% dataMod = qammod(input_data_symbols, M, 'bin'); % modulate the input data using binary encoding
dataMod_GrayCode256 = qammod(input_data_symbols256, M256, 'UnitAveragePower', false); % ... and Gray coded
avg_mod_pwr256 = sum(abs(dataMod_GrayCode256).^2)/M256

EbN0256 = 10;
snr0_256 = EbN0256 + 10*log10(k256) - 10*log10(samplesPerSymbol256);
receivedSignal_GC256 = awgn(dataMod_GrayCode256, snr0_256, 'measured');
s_plot = scatterplot(receivedSignal_GC256, 1, 0, 'g.');
hold on
scatterplot(dataMod_GrayCode256, 1, 0, 'k*', s_plot);

for i = 1:50
    snr256(i) = EbN0256 + 10*log10(k256) - 10*log10(samplesPerSymbol256);
    receivedSignal_GrayCode256 = awgn(dataMod_GrayCode256, snr256(i), 'measured');
    avg_receivedPwr256 = sum(abs(receivedSignal_GrayCode256).^2)/M256
    output_data_symbols_GrayCode256 = qamdemod(receivedSignal_GrayCode256, M256);

    output_data_matrix_GrayCode256 = de2bi(output_data_symbols_GrayCode256, k256);
    output_data_GrayCode256 = output_data_matrix_GrayCode256(:);

    [numbErrors_GrayCode256, BER_GrayCode256] = biterr(input_data256, output_data_GrayCode256);
    logber256(i) = log2(BER_GrayCode256);
    ber256(i) = BER_GrayCode256;
    samplesPerSymbol256 = samplesPerSymbol256 + 1;
end

figure, plot(snr256, logber256); xlabel('SNR'), ylabel('log(BER)'), title('log(BER) vs SNR, 256-QAM (adjusted samples per symbol');

