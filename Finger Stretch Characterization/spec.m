function spec(x, Fs)
    window = hamming(512);
    noverlap = 256;
    nfft = 4096;
    spectrogram(x, window, noverlap, nfft, Fs);
    ylabel('Time (s)');
    xlabel('Frequency kHz');
    colormap("parula");
    colorbar;
end