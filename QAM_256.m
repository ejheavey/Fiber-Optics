M = 256;
k = log2(M);
n = 30000;
samplesPerSymbol = 1;

rng default;

input_data = randi([0 1], n, 1);
input_data_matrix = reshape(input_data, length(input_data)/k, k); % reshape the data to n x n 
input_data_symbols = bi2de(input_data_matrix); % binary to decimal conversion
figure; stem(input_data_symbols(1:10)); % show the input symbols in decimal

% dataMod = qammod(input_data_symbols, M, 'bin'); % modulate the input data using binary encoding
dataMod_GrayCode = qammod(input_data_symbols, M, 'UnitAveragePower', false); % ... and Gray coded
avg_mod_pwr = sum(abs(dataMod_GrayCode).^2)/M

EbN0 = 10;
snr = EbN0 + 10*log10(k) - 10*log10(samplesPerSymbol) % EbN0 + signalPower - noisePower

% receivedSignal = awgn(dataMod, snr, 'measured');
receivedSignal_GrayCode = awgn(dataMod_GrayCode, snr, 'measured');
avg_receivedPwr = sum(abs(receivedSignal_GrayCode).^2)/M

s_plot = scatterplot(receivedSignal, 1, 0, 'g.');
hold on
scatterplot(dataMod, 1, 0, 'k*', s_plot);

output_data_symbols = qamdemod(receivedSignal, M, 'bin');
output_data_symbols_GrayCode = qamdemod(receivedSignal_GrayCode, M);

output_data_matrix = de2bi(output_data_symbols, k);
output_data = output_data_matrix(:);

output_data_matrix_GrayCode = de2bi(output_data_symbols_GrayCode, k);
output_data_GrayCode = output_data_matrix_GrayCode(:);

[numbErrors, BER] = biterr(input_data, output_data);
[numbErrors_GrayCode, BER_GrayCode] = biterr(input_data, output_data_GrayCode);
fprintf('\nBinary Coding BER: %5.2e, based on %d errors\n', BER, numbErrors);
fprintf('\nGray Coded BER: %5.2e, based on %d errors\n', BER_GrayCode, numbErrors_GrayCode);

x = (0:255); % Integer input
y1 = qammod(x, 256, 'bin');

scatterplot(y1);
text(real(y1)+0.1, imag(y1), dec2bin(x));
title('16-QAM, Binary Symbol Mapping');
axis([-4 4 -4 4])

%% Gray-coded symbol mapping
y2 = qammod(x, 256, 'gray');

scatterplot(y2)
text(real(y2)+0.1, imag(y2), dec2bin(x))
title('256-QAM, Gray-coded Symbol Mapping')
axis([-16 16 -16 16])



