Fs = 96;  % Sampling Frequency
LFpass = 3; 
LFstop = 5; 
LDpass = 0.057501127785;  
LDstop = 0.0001; 
Ldens  = 20;   
[N, Fo, Ao, W] = firpmord([LFpass, LFstop]/(Fs/2), [1 0], [LDpass, LDstop]); 
Lowpass  = firpm(N, Fo, Ao, W, {Ldens});

HFstop = 2;  
HFpass = 4; 
HAstop = 80;  
HApass = 1; 
Hmatch = 'stopband';  
h  = fdesign.highpass(HFstop, HFpass, HAstop, HApass, Fs);
Highpass = design(h, 'butter', 'MatchExactly', Hmatch);
     
% AM Modified
figure(1);subplot(2,1,1);plot(Lowpass);title('LowPass');
figure(1);subplot(2,1,2);plot(Highpass);title('HighPass');

rec = 20;
% AM Modified
% [time, X, Y, Z, P, R] = openrecallxyz(rec,0, Start_time, End_time);
% plot(time, Y)
% xlabel('time (s)'); ylabel('normalised force'); 
% [time, X, Y, Z, P, R] = openrecallxyz(rec,0, Start_time, End_time);
rawDataPath = 'C:\Users\engs1602\research\data\Data Feb-Mar 2016\Raw Data - wiimote';
file = fopen(fullfile(rawDataPath,sprintf('%d.txt',rec)),'rt');
raw = textscan(file, '%f:%f:%f:%f \t %f \t %f \t %f \t %f \t %f \t');
fclose(file);
time = raw{4}/1000 + raw{3} + 60*raw{2};
Y = raw{6}; 
Z = raw{7};
Y = Y - mean(Y);
Z = Z - mean(Z);

    % Filter signals
    cliplength = 256; overlap = 128;
    
    Y_s = filtfilt(Lowpass, 1, Y);
    Z_s = filtfilt(Lowpass,1,Z-mean(Z));
    arc = atand(Y_s./(Z_s+mean(Z)));
    
    % AM Modified
    figure(rec);subplot(4,1,1);plot(time,Y);ylabel('Y');
    figure(rec);subplot(4,1,2);plot(time,Y_s);ylabel('Y Low-passed');
        
    % AM Modified
    % Since sortoutpeaks3.m N/A, interval2 hard-coded by eye-balling for now
    % Just removes stationary periods at start and end of recording    
    % [interval2, peak_locs, trough_locs] = sortoutpeaks3(Y_s, 0.1);
%     [interval2] = sortoutpeaks3(Y_s, 0.1);   
    interval2 = sortoutpeaks3FarahNA(Y_s, 0.1); 
    time = time(interval2);
    Y = Y(interval2);
    % AM Modified
    figure(rec);subplot(4,1,3);plot(time,Y);ylabel('Y ends removed');
   
    idx = 1;
    intstart = 1;
    while intstart < length(interval2) - cliplength
        intstart;
        intend = intstart + cliplength;
        % AM Comments
        % The following line does not make sense!?
        Yf = [Yf, filter(Highpass, Y(intstart:intend))];%***
        Yfhat = fft(Yf, cliplength);
        spectra20(:,idx) = abs(Yfhat(1:cliplength/2));
        %rec_no(idx) = rec;
        spd(idx) = sum(abs(diff(arc(intstart:intend))));
        idx = idx+1;
        intstart = intstart + overlap;
    end

    figure();
    subplot(3,1,1:2)
    plot(spectra20(:, 10:15), 'b')
    hold on
    plot(median(spectra20(:, :), 2), 'LineWidth', 3)
    hold off
    plot(spectra(:, 10:15), 'b')   
    xlabel('Frequency'); ylabel('Magnitude')
   
subtightplot(3,1,1:2)
imagesc(spectra20(10:end, :)); caxis([0, 2.5])
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'Fontsize', 12);
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'Fontsize', 12);
set(gca,'YDir','normal'); ylabel('frequency')
subtightplot(3,1,3)
plot(spd); xlim([1,73]); ylim([150, 230])
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'Fontsize', 12);
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'Fontsize', 12);
xlabel('time'); ylabel('speed')


subtightplot(2,1,1)
plot(time(interval2), Y); xlim([37, 44])
subtightplot(2,1,2)
plot(time(interval2),  filter(Highpass, Y)); xlim([37, 44])
% Collect all spectra for pump 1

freqs = [14,21,46,60];
for f = 1:3
    a = spectra114(freqs(f), :);
    a = a(a~=0)';
    figure(10);
    subplot(4,1,f)
    histfit(a,40,'gamma'); xlim([0, 4]); ylim([0, 300])
%     set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'Fontsize', 10);
    title(sprintf('Frequency = %d', freqs(f)))
end
f = 4;
a = spectra114(freqs(f), :);
a = a(a~=0)';
subplot(4,1,f)
histfit(a,40,'gamma'); xlim([0, 4]); ylim([0, 300])
set(gca, 'Fontsize', 10);
title(sprintf('Frequency = %d', freqs(f)))

fac = [2,3,6,7];
for f = 1:64
    a = spectra114(f, :);
    a = a(a~=0)';
    % AM Modified
%     subtightplot(8,8,f, [0.02, 0.01])
    figure(11);
    subplot(8,8,f);
    histfit(a,100,'gamma'); xlim([0,8])
end

for f = 1:4
    a = spectra114(fac(f)*8-2, :);
    a = a(a~=0)';
    figure(12);
    subplot(4,1,f)
    histfit(a,100,'gamma'); xlim([0,5])
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'Fontsize', 10);
    title(sprintf('Frequency = %d', fac(f)*8-2))
end
xlabel('Magnitude')
idx = 1;

for rec = 1:16

file = fopen(strcat(num2str(rec), '.txt'), 'rt');
raw = textscan(file, '%f:%f:%f:%f \t %f \t %f \t %f \t %f \t %f \t');
fclose(file);
time = raw{4}/1000 + raw{3} + 60*raw{2};
Y = raw{6}; Z = raw{7};
Y = Y - mean(Y);
plot(time, Y)

    Y_s = filtfilt(Lowpass, 1, Y);
    Z_s = filtfilt(Lowpass,1,Z-mean(Z));
    arc = atand(Y_s./(Z_s+mean(Z)));
    %     [interval2, peak_locs, trough_locs] = sortoutpeaks3(Y_s, 0.1);
    [interval2] = sortoutpeaks3(Y_s, 0.1);    Y = Y(interval2);
    
    intstart = 1;
    while intstart < length(interval2) - cliplength
        intstart;
        intend = intstart + cliplength;
        Yf = filter(Highpass, Y(intstart:intend));
        Yfhat = fft(Yf, cliplength);
        spectra(:,idx) = abs(Yfhat(1:cliplength/2));
        rec_no(idx) = rec;
        spd(idx) = sum(abs(diff(arc(intstart:intend))));
        idx = idx+1;
        intstart = intend;
    end
end

figure()
imagesc(log(spectra(6:end, :))); caxis([-4, 1])
hist(log(spectra(:)), 100)

for f = 1:64
    a = spectra(f, :);
    a = a(a~=0)';
    subtightplot(8,8,f, [0.02, 0.01])
    histfit(a,100,'gamma'); 
end

HFstop = 3;  HFpass = 6; HAstop = 80;  HApass = 1; Hmatch = 'stopband';  
h  = fdesign.highpass(HFstop, HFpass, HAstop, HApass, Fs);
Highpass = design(h, 'butter', 'MatchExactly', Hmatch);

more = []
less = []

rec = 32; idx = 1;
file = fopen(strcat(num2str(rec), '.txt'), 'rt');
raw = textscan(file, '%f:%f:%f:%f \t %f \t %f \t %f \t %f \t %f \t');
fclose(file);
time = raw{4}/1000 + raw{3} + 60*raw{2};
Y = raw{6}; Z = raw{7};
Y = Y - mean(Y);
plot(time, Y)

   Y_s = filtfilt(Lowpass, 1, Y);
    Z_s = filtfilt(Lowpass,1,Z-mean(Z));
    arc = atand(Y_s./(Z_s+mean(Z)));
    %     [interval2, peak_locs, trough_locs] = sortoutpeaks3(Y_s, 0.1);
    [interval2] = sortoutpeaks3(Y_s, 0.1);    Y = Y(interval2);
    
    intstart = 1; 
    while intstart < length(interval2) - cliplength
        intstart;
        intend = intstart + cliplength;
        Yf = filter(Highpass, Y(intstart:intend));
        %Yfhat = fft(Y(intstart:intend), cliplength); %****
        Yfhat = fft(Yf, cliplength); 
        spectra27(:,idx) = abs(Yfhat(1:cliplength/2));
        %rec_no(idx) = rec;
        spd27(idx) = sum(abs(diff(arc(intstart:intend))));
        idx = idx+1;
        intstart = intstart+overlap;
    end


plot(median(spectra27(:, 1:62), 2), 'g')
hold on
plot(median(spectra27(:, 63:end), 2), 'b')
hold off

specg = [zeros(62, 1); ones(122, 1)];
for f = 1:64
    boxplot(spectra27(f, :), specg)
    pause(0.5)
end

 Yf = filter(Highpass, Y);
    plot(time, [Y, Yf])
imagesc(log(spectra27))

idx = 1;
for rec = 66:84 
file = fopen(strcat(num2str(rec), '.txt'), 'rt');
raw = textscan(file, '%f:%f:%f:%f \t %f \t %f \t %f \t %f \t %f \t');
fclose(file);
time = raw{4}/1000 + raw{3} + 60*raw{2};
Y = raw{6}; Z = raw{7};
Y = Y - mean(Y);
plot(time, Y)

   Y_s = filtfilt(Lowpass, 1, Y);
    Z_s = filtfilt(Lowpass,1,Z-mean(Z));
    arc = atand(Y_s./(Z_s+mean(Z)));
    %     [interval2, peak_locs, trough_locs] = sortoutpeaks3(Y_s, 0.1);
    [interval2] = sortoutpeaks3(Y_s, 0.1);    Y = Y(interval2);
    
    intstart = 1; 
    while intstart < length(interval2) - cliplength
        intstart;
        intend = intstart + cliplength;
        Yf = filter(Highpass, Y(intstart:intend));
        %Yfhat = fft(Y(intstart:intend), cliplength); %****
        Yfhat = fft(Yf, cliplength); 
        spectra18(:,idx) = abs(Yfhat(1:cliplength/2));
        rec_no18(idx) = rec;
        spd18(idx) = sum(abs(diff(arc(intstart:intend))));
        idx = idx+1;
        intstart = intstart+overlap;
    end
end

imagesc(spectra18); caxis([0,4])
plot(median(spectra18(:, rec_no18 <=71), 2), 'g', 'LineWidth', 6)
hold on
plot(median(spectra18(:, rec_no18 >= 72), 2), 'b', 'LineWidth', 6)
hold off

hold on
for r = 66:71
    plot(median(spectra18(:, rec_no18 == r), 2), 'g--')
end
for r = 72:84
    plot(median(spectra18(:, rec_no18 == r), 2), 'b--')
end
hold off
    plot(median(spectra18(:, rec_no18 <=71), 2), 'g')
    
xlabel('Frequency')
ylabel('Magnitude')
    

for f = 1:4
    a = spectra(f*16, rec_no18 <=71);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h, x] = hist(a, 30);
      g = pdf(pd, x);
    subplot(4,1,f)
      plot(x, g, 'g')
      
        a = spectra(f*16, rec_no18 >=72);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h, x] = hist(a, 30);
      g = pdf(pd, x);
    %subtightplot(8,8,f, [0.02, 0.01])
    hold on
      plot(x, g, 'b'); xlim([0, 1.5])
      hold off
      set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'Fontsize', 12);
    ylabel('PDF')
    title(sprintf('Frequency = %d', f*16))
end
xlabel('Magnitude')

for f = 1:64
    a = spectra(f, rec_no18 <=71);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h1, x1] = hist(a, 30);
      g1 = pdf(pd, x1);
    subtightplot(8,8,f, [0.02, 0.01])
      bar(x1, g1, 'g')
      
        a = spectra(f, rec_no18 >=72);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h2, x2] = hist(a, 30);
      g2 = pdf(pd, x2);
    subtightplot(8,8,f, [0.02, 0.01])
    hold on
      plot(x2, g2, 'b')
      hold off
end

%% Jabalini then vs now

%Now
idx = 1
a = 1:85; 
for rec = 303:313;
    rec
    file = fopen(strcat(num2str(rec), '.txt'), 'rt');
raw = textscan(file, '%f:%f:%f:%f \t %f \t %f \t %f \t %f \t %f \t');
fclose(file);
time = raw{4}/1000 + raw{3} + 60*raw{2};
Y = raw{6}; Z = raw{7};
Y = Y - mean(Y);
    %
    %     all_peak_locs=[];
    %     all_trough_locs=[];
    %
    %     %[interval2, peak_locs, trough_locs] = sortoutpeaks3(Ey, 0.1);
    Y_s = filtfilt(Lowpass, 1, Y);
    Z_s = filtfilt(Lowpass,1,Z-mean(Z));
    arc = atand(Y_s./(Z_s+mean(Z)));
    %     [interval2, peak_locs, trough_locs] = sortoutpeaks3(Y_s, 0.1);
    [interval2] = sortoutpeaks3(Y_s, 0.1);    Y = Y(interval2);
    
    intstart = 1;
    while intstart < length(interval2) - cliplength
        intend = intstart + cliplength;
        %Ys(:, idx) = Y(intstart:intend);
        Yf = filter(Highpass, Y(intstart:intend)); 
        Yfhat = fft(Yf, cliplength);
        spectra(:,idx) = abs(Yfhat(1:cliplength/2));
        rec_no(idx) = rec;
        spd(idx) = sum(abs(diff(arc(intstart:intend))));
        idx = idx+1;
        intstart = intstart+overlap;
    end
end

spectra16 = spectra;

freqs = [10,22,42,62];
 for f = 1:4
    a = spectra16(freqs(f), :);
    a = a(a~=0);
    pd = fitdist(a', 'Gamma');
    [h, x] = hist(a, 50);
    g = pdf(pd, x);
    %subtightplot(6, 4, (f-1)/5+1)
    subplot(2,2,f)
    %bar(x, h/sum(h)/(x(2)-x(1))); %xlim([0,4])
    hold on
    plot(x, g, 'b')
    %hold off
    
    a = spectra14(freqs(f), :);
    a = a(a~=0);
    pd = fitdist(a', 'Gamma');
    [h, x] = hist(a, 50);
    g = pdf(pd, x);
    
    %bar(x, h/sum(h)/(x(2)-x(1))); %xlim([0,4])
    %hold on
    plot(x, g, 'g')
    hold off
    xlim([0,5])
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'Fontsize', 12);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'Fontsize', 12);
    ylabel('PDF');  xlabel('Magnitude')
    title(sprintf('Frequency = %d', freqs(f)))
 end
 xlabel('Magnitude')
 
imagesc(spectra18); caxis([0,4])
plot((48/128)*(1:128), median(spectra14(:, :), 2), 'g', 'LineWidth', 2)
hold on
plot((48/128)*(1:128), prctile(spectra14(:, :)', 95), 'g', 'LineWidth', 1)
plot((48/128)*(1:128), median(spectra16(:, :), 2), 'b', 'LineWidth', 2)
plot((48/128)*(1:128), prctile(spectra16(:, :)', 95), 'b', 'LineWidth', 1)
hold off
xlim([10, 48])
xlabel('Frequency', 'fontsize',16); ylabel('Magnitude', 'fontsize',16);

hold on
for r = 66:71
    plot(median(spectra14(:, rec_no18 == r), 2), 'g--')
end
for r = 72:84
    plot(median(spectra16(:, rec_no18 == r), 2), 'b--')
end
hold off
    plot(median(spectra18(:, rec_no18 <=71), 2), 'g')

    plot(spectra18(:, 100:102))
    
    
 %% From other script
 
 idx = 1;

for rec = 1:16

file = fopen(strcat(num2str(rec), '.txt'), 'rt');
raw = textscan(file, '%f:%f:%f:%f \t %f \t %f \t %f \t %f \t %f \t');
fclose(file);
time = raw{4}/1000 + raw{3} + 60*raw{2};
Y = raw{6}; Z = raw{7};
Y = Y - mean(Y);
plot(time, Y)

    Y_s = filtfilt(Lowpass, 1, Y);
    Z_s = filtfilt(Lowpass,1,Z-mean(Z));
    arc = atand(Y_s./(Z_s+mean(Z)));
    %     [interval2, peak_locs, trough_locs] = sortoutpeaks3(Y_s, 0.1);
    [interval2] = sortoutpeaks3(Y_s, 0.1);    Y = Y(interval2);
    
    intstart = 1;
    while intstart < length(interval2) - cliplength
        intstart;
        intend = intstart + cliplength;
        Yf = filter(Highpass, Y(intstart:intend));
        Yfhat = fft(Yf, cliplength);
        spectra(:,idx) = abs(Yfhat(1:cliplength/2));
        rec_no(idx) = rec;
        spd(idx) = sum(abs(diff(arc(intstart:intend))));
        idx = idx+1;
        intstart = intend;
    end
end

imagesc(log(spectra(6:end, :))); caxis([-4, 1])
hist(log(spectra(:)), 100)

for f = 1:64
    a = spectra(f, :);
    a = a(a~=0)';
    subtightplot(8,8,f, [0.02, 0.01])
    histfit(a,100,'gamma'); 
end

HFstop = 3;  HFpass = 6; HAstop = 80;  HApass = 1; Hmatch = 'stopband';  
h  = fdesign.highpass(HFstop, HFpass, HAstop, HApass, Fs);
Highpass = design(h, 'butter', 'MatchExactly', Hmatch);

more = []
less = []

rec = 32; idx = 1;
file = fopen(strcat(num2str(rec), '.txt'), 'rt');
raw = textscan(file, '%f:%f:%f:%f \t %f \t %f \t %f \t %f \t %f \t');
fclose(file);
time = raw{4}/1000 + raw{3} + 60*raw{2};
Y = raw{6}; Z = raw{7};
Y = Y - mean(Y);
plot(time, Y)

    Y_s = filtfilt(Lowpass, 1, Y);
    Z_s = filtfilt(Lowpass,1,Z-mean(Z));
    arc = atand(Y_s./(Z_s+mean(Z)));
    %     [interval2, peak_locs, trough_locs] = sortoutpeaks3(Y_s, 0.1);
    [interval2] = sortoutpeaks3(Y_s, 0.1);    Y = Y(interval2);
    
    intstart = 1; 
    while intstart < length(interval2) - cliplength
        intstart;
        intend = intstart + cliplength;
        Yf = filter(Highpass, Y(intstart:intend));
        %Yfhat = fft(Y(intstart:intend), cliplength); %****
        Yfhat = fft(Yf, cliplength); 
        spectra27(:,idx) = abs(Yfhat(1:cliplength/2));
        %rec_no(idx) = rec;
        spd27(idx) = sum(abs(diff(arc(intstart:intend))));
        idx = idx+1;
        intstart = intstart+overlap;
    end


plot(median(spectra27(:, 1:62), 2), 'g')
hold on
plot(median(spectra27(:, 63:end), 2), 'b')
hold off

specg = [zeros(62, 1); ones(122, 1)];
for f = 1:64
    boxplot(spectra27(f, :), specg)
    pause(0.5)
end

 Yf = filter(Highpass, Y);
    plot(time, [Y, Yf])
imagesc(log(spectra27))

idx = 1;
for rec = 66:84 
file = fopen(strcat(num2str(rec), '.txt'), 'rt');
raw = textscan(file, '%f:%f:%f:%f \t %f \t %f \t %f \t %f \t %f \t');
fclose(file);
time = raw{4}/1000 + raw{3} + 60*raw{2};
Y = raw{6}; Z = raw{7};
Y = Y - mean(Y);
plot(time, Y)

   Y_s = filtfilt(Lowpass, 1, Y);
    Z_s = filtfilt(Lowpass,1,Z-mean(Z));
    arc = atand(Y_s./(Z_s+mean(Z)));
    %     [interval2, peak_locs, trough_locs] = sortoutpeaks3(Y_s, 0.1);
    [interval2] = sortoutpeaks3(Y_s, 0.1);    Y = Y(interval2);
    
    intstart = 1; 
    while intstart < length(interval2) - cliplength
        intstart;
        intend = intstart + cliplength;
        Yf = filter(Highpass, Y(intstart:intend));
        %Yfhat = fft(Y(intstart:intend), cliplength); %****
        Yfhat = fft(Yf, cliplength); 
        spectra18(:,idx) = abs(Yfhat(1:cliplength/2));
        rec_no18(idx) = rec;
        spd18(idx) = sum(abs(diff(arc(intstart:intend))));
        idx = idx+1;
        intstart = intstart+overlap;
    end
end

imagesc(spectra18); caxis([0,4])
plot(median(spectra18(:, rec_no18 <=71), 2), 'g', 'LineWidth', 6)
hold on
plot(median(spectra18(:, rec_no18 >= 72), 2), 'b', 'LineWidth', 6)
hold off

hold on
for r = 66:71
    plot(median(spectra18(:, rec_no18 == r), 2), 'g--')
end
for r = 72:84
    plot(median(spectra18(:, rec_no18 == r), 2), 'b--')
end
hold off
    plot(median(spectra18(:, rec_no18 <=71), 2), 'g')
    
xlabel('Frequency')
ylabel('Magnitude')
    

for f = 1:4
    a = spectra(f*16, rec_no18 <=71);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h, x] = hist(a, 30);
      g = pdf(pd, x);
    subplot(4,1,f)
      plot(x, g, 'g')
      
        a = spectra(f*16, rec_no18 >=72);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h, x] = hist(a, 30);
      g = pdf(pd, x);
    %subtightplot(8,8,f, [0.02, 0.01])
    hold on
      plot(x, g, 'b'); xlim([0, 1.5])
      hold off
      set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'Fontsize', 12);
    ylabel('PDF')
    title(sprintf('Frequency = %d', f*16))
end
xlabel('Magnitude')

for f = 1:64
    a = spectra(f, rec_no18 <=71);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h1, x1] = hist(a, 30);
      g1 = pdf(pd, x1);
    subtightplot(8,8,f, [0.02, 0.01])
      bar(x1, g1, 'g')
      
        a = spectra(f, rec_no18 >=72);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h2, x2] = hist(a, 30);
      g2 = pdf(pd, x2);
    subtightplot(8,8,f, [0.02, 0.01])
    hold on
      plot(x2, g2, 'b')
      hold off
end

for f = 6:5:64
    for  loc = 1:11
        subtightplot(11,1,loc)
            hist(spectra(f, loc_num(rec_no) == loc),20); 
            thresh = prctile(spectra(f, loc_num(rec_no) == loc), 95);
            max = get(gca, 'YLim');
            hold on
            plot([thresh, thresh], [0,max(2)*0.8], 'r')
            hold off
            if loc ~= 11
                set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'Fontsize', 12);
            else
                xlabel('Magnitude')
            end
            set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'Fontsize', 12);
            xlim([0,4])
    end
    suptitle(sprintf('%d%', f)) 
    pause(1)
end

%% Paper

lims = prctile(lspectra(:), [5,95]); 
for l = 1:29 subplot(3,1,1:2) 
    imagesc(lspectra(:, loc_num(rec_no) == l)) 
    caxis(lims); set(gca,'YDir','normal') x = get(gca, 'XLim'); 
    title(['Pump', char(32), num2str(l)]) subplot(3,1,3) plot(spd(loc_num(rec_no) == l)); 
    xlim(x) 
    print -dpsc '-append' 'ken16spectraspeed' 
end

l = 3;
imagesc(lspectra(10:end, loc_num(rec_no) == l)) ; caxis(prctile(reshape(lspectra(10:end, loc_num(rec_no) == l), 1, []), [5,95]));set(gca,'YDir','normal')
recdiff = diff(rec_no(loc_num(rec_no) == l)); ends = find(recdiff ~= 0);
hold on
for i = 1:length(ends)
    plot([ends(i),ends(i)], [0, 128], 'k', 'LineWidth', 2)
end
hold off
    
%% doing histograms


l = 19;

spec19 = spectra(1:end, loc_num(rec_no) == l);
rec19 = rec_no(loc_num(rec_no) == l);
for fi = 1:64
    f = 2*fi;
    a = spec19(f, rec19 <=190);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h1, x1] = hist(a, 30);
      g1 = pdf(pd, x1);
    subtightplot(8,8,fi, [0.02, 0.01])
      bar(x1, h1, 'g')
      
        a = spec19(f, rec19 >=191);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h2, x2] = hist(a, 30);
      g2 = pdf(pd, x2);
    subtightplot(8,8,fi, [0.02, 0.01])
    hold on
      bar(x2, h2, 'b')
      hold off
end

plot((48/128)*(1:128), median(spec19(:,rec19 <=190), 2), 'g','LineWidth', 2);
hold on
plot((48/128)*(1:128), prctile(spec19(:,rec19 <=190)', [95])', 'g', 'LineWidth', 1)
plot((48/128)*(1:128), median(spec19(:,rec19 >190), 2), 'b', 'LineWidth', 2)
plot((48/128)*(1:128), prctile(spec19(:,rec19 >190)', [95])', 'b', 'LineWidth', 1)
hold off
xlim([10, 48])
xlabel('Frequency'); ylabel('Magnitude');


l = 19;

spec19 = spectra(1:end, loc_num16(rec_no) == l);
rec19 = rec_no(loc_num16(rec_no) == l);
for fi = 1:64
    f = 2*fi;
    a = spec19(f, rec19 <=190);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h1, x1] = hist(a, 30);
      g1 = pdf(pd, x1);
    subtightplot(8,8,fi, [0.02, 0.01])
      bar(x1, h1, 'g')
      
        a = spec19(f, rec19 >=191);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h2, x2] = hist(a, 30);
      g2 = pdf(pd, x2);
    subtightplot(8,8,fi, [0.02, 0.01])
    hold on
      bar(x2, h2, 'b')
      hold off
end

plot((48/128)*(1:128), median(spec19(:,rec19 <=190), 2), 'g','LineWidth', 2);
hold on
plot((48/128)*(1:128), prctile(spec19(:,rec19 <=190)', [95])', 'g', 'LineWidth', 1)
plot((48/128)*(1:128), median(spec19(:,rec19 >190), 2), 'b', 'LineWidth', 2)
plot((48/128)*(1:128), prctile(spec19(:,rec19 >190)', [95])', 'b', 'LineWidth', 1)
hold off
xlim([10, 48])
xlabel('Frequency', 'fontsize',16); ylabel('Magnitude', 'fontsize',16);
'fontsize',20

l = 3;

spec3 = spectra(1:end, loc_num(rec_no) == l);
imagesc(log(spec3)); caxis(prctile(log(spec3(:)), [5,95]))
rec3 = rec_no(loc_num(rec_no) == l);
for fi = 1:64
    f = 2*fi;
    a = spec3(f, rec3 <=71);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h1, x1] = hist(a, 30);
      g1 = pdf(pd, x1);
    subtightplot(8,8,fi, [0.02, 0.01])
      bar(x1, h1, 'g')
      
        a = spec3(f, rec3 >=191);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h2, x2] = hist(a, 30);
      g2 = pdf(pd, x2);
    subtightplot(8,8,fi, [0.02, 0.01])
    hold on
      bar(x2, h2, 'b')
      hold off
end

plot((48/128)*(1:128), median(spec3(:,rec19 <=190), 2), 'g','LineWidth', 2);
hold on
plot((48/128)*(1:128), prctile(spec3(:,rec19 <=190)', [95])', 'g', 'LineWidth', 1)
plot((48/128)*(1:128), median(spec3(:,rec19 >190), 2), 'b', 'LineWidth', 2)
plot((48/128)*(1:128), prctile(spec3(:,rec19 >190)', [95])', 'b', 'LineWidth', 1)
hold off
xlim([10, 48])
xlabel('Frequency'); ylabel('Magnitude');

l = 6;

spec6 = spectra(1:end, loc_num16(rec_no) == l);
rec6 = rec_no(loc_num16(rec_no) == l);
for fi = 1:64
    f = 2*fi;
    a = spec6(f, rec6 <=71);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h1, x1] = hist(a, 30);
      g1 = pdf(pd, x1);
    subtightplot(8,8,fi, [0.02, 0.01])
      plot(x1, h1/sum(rec6 <=71), 'g')
      
        a = spec6(f, rec6 >=72);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h2, x2] = hist(a, 30);
      g2 = pdf(pd, x2);
    subtightplot(8,8,fi, [0.02, 0.01])
    hold on
      plot(x2, h2/sum(rec6 >=72), 'b')
      hold off
end

plot((48/128)*(1:128), median(spec6(:,rec6 <=71), 2), 'g','LineWidth', 2);
hold on
plot((48/128)*(1:128), prctile(spec6(:,rec6 <=71)', [95])', 'g', 'LineWidth', 1)
plot((48/128)*(1:128), median(spec6(:,rec6 >71), 2), 'b', 'LineWidth', 2)
plot((48/128)*(1:128), prctile(spec6(:,rec6 >71)', [95])', 'b', 'LineWidth', 1)
hold off
xlim([10, 48])
xlabel('Frequency', 'fontsize',16); ylabel('Magnitude', 'fontsize',16);

%% Do Jab then/now again
idx = 1;
Fs = 96;  % Sampling Frequency


LFpass = 3; LFstop = 5; LDpass = 0.057501127785;  LDstop = 0.0001; Ldens  = 20;   
[N, Fo, Ao, W] = firpmord([LFpass, LFstop]/(Fs/2), [1 0], [LDpass, LDstop]); Lowpass  = firpm(N, Fo, Ao, W, {Ldens});

HFstop = 2;  HFpass = 4; HAstop = 80;  HApass = 1; Hmatch = 'stopband';  
h  = fdesign.highpass(HFstop, HFpass, HAstop, HApass, Fs);
Highpass = design(h, 'butter', 'MatchExactly', Hmatch);
idx = 1;
cliplength = 256; overlap = 128;



idx = 1;
a = 1:85; 
for rec = a(loc_num == 1);
    [time, X, Y, Z] = openrecallxyz(rec,0, Start_time, End_time);
    %
    %     all_peak_locs=[];
    %     all_trough_locs=[];
    %
    %     %[interval2, peak_locs, trough_locs] = sortoutpeaks3(Ey, 0.1);
    Y_s = filtfilt(Lowpass, 1, Y);
    Z_s = filtfilt(Lowpass,1,Z-mean(Z));
    arc = atand(Y_s./(Z_s+mean(Z)));
    %     [interval2, peak_locs, trough_locs] = sortoutpeaks3(Y_s, 0.1);
    [interval2] = sortoutpeaks3(Y_s, 0.1);    Y = Y(interval2);
    
    intstart = 1;
     Yf = filter(Highpass, Y);
    while intstart < length(interval2) - cliplength
        intend = intstart + cliplength;
        %Ys(:, idx) = Y(intstart:intend);
       Yfhat = fft(Yf(intstart:intend), cliplength);
        spectra(:,idx) = abs(Yfhat(1:cliplength/2));
        rec_no(idx) = rec;
        spd(idx) = sum(abs(diff(arc(intstart:intend))));
        idx = idx+1;
        intstart = intstart+overlap;
    end
end

idx = 1;
for rec = 303:313;
    rec
    file = fopen(strcat(num2str(rec), '.txt'), 'rt');
raw = textscan(file, '%f:%f:%f:%f \t %f \t %f \t %f \t %f \t %f \t');
fclose(file);
time = raw{4}/1000 + raw{3} + 60*raw{2};
Y = raw{6}; Z = raw{7};
Y = Y - mean(Y);
    %
    %     all_peak_locs=[];
    %     all_trough_locs=[];
    %
    %     %[interval2, peak_locs, trough_locs] = sortoutpeaks3(Ey, 0.1);
    Y_s = filtfilt(Lowpass, 1, Y);
    Z_s = filtfilt(Lowpass,1,Z-mean(Z));
    arc = atand(Y_s./(Z_s+mean(Z)));
    %     [interval2, peak_locs, trough_locs] = sortoutpeaks3(Y_s, 0.1);
    [interval2] = sortoutpeaks3(Y_s, 0.1);    Y = Y(interval2);
    
    intstart = 1;
    while intstart < length(interval2) - cliplength
        intend = intstart + cliplength;
        %Ys(:, idx) = Y(intstart:intend);
        Yf = filter(Highpass, Y(intstart:intend)); 
        Yfhat = fft(Yf, cliplength);
        spectra(:,idx) = abs(Yfhat(1:cliplength/2));
        rec_no(idx) = rec;
        spd(idx) = sum(abs(diff(arc(intstart:intend))));
        idx = idx+1;
        intstart = intstart+overlap;
    end
end

spectra16 = spectra;

% Final plot
fs = 10; fs1 = 12;


subplot(3,1,1)
plot((48/128)*(1:128), median(spec19(:,rec19 <=190), 2), 'b','LineWidth', 2);
hold on
plot((48/128)*(1:128), prctile(spec19(:,rec19 <=190)', [95])', 'b', 'LineWidth', 1)
plot((48/128)*(1:128), median(spec19(:,rec19 >190), 2), 'g', 'LineWidth', 2)
plot((48/128)*(1:128), prctile(spec19(:,rec19 >190)', [95])', 'g', 'LineWidth', 1)
hold off
xlim([10, 48])
xlabel('Frequency (Hz)', 'fontsize',fs); ylabel('Magnitude', 'fontsize',fs);
title('(a)', 'fontsize',fs1)
subplot(3,1,2)


plot((48/128)*(1:128), median(spec6(:,rec6 <=71), 2), 'b','LineWidth', 2);
hold on
plot((48/128)*(1:128), prctile(spec6(:,rec6 <=71)', [95])', 'b', 'LineWidth', 1)
plot((48/128)*(1:128), median(spec6(:,rec6 >71), 2), 'g', 'LineWidth', 2)
plot((48/128)*(1:128), prctile(spec6(:,rec6 >71)', [95])', 'g', 'LineWidth', 1)
hold off
xlim([10, 48])
xlabel('Frequency (Hz)', 'fontsize',fs); ylabel('Magnitude', 'fontsize',fs);
title('(b)', 'fontsize',fs1)
subplot(3,1,3)
plot((48/128)*(1:128), median(spectra14(:, :), 2), 'g', 'LineWidth', 2)
hold on
plot((48/128)*(1:128), prctile(spectra14(:, :)', 95), 'g', 'LineWidth', 1)
plot((48/128)*(1:128), median(spectra16(:, :), 2), 'b', 'LineWidth', 2)
plot((48/128)*(1:128), prctile(spectra16(:, :)', 95), 'b', 'LineWidth', 1)
hold off
xlim([10, 48])
xlabel('Frequency (Hz)', 'fontsize',fs); ylabel('Magnitude', 'fontsize',fs);
title('(c)', 'fontsize',fs1)
int = 0.3
% histogram
[h1, x1 ] =  hist(spec19(61,rec19 <=190), 15);
bh1 = bar(x1, h1/sum(h1)/(x1(2)-x1(1)), 1)

set(bh1,'facecolor','b');
set(get(bh1,'Children'),'FaceAlpha',0.2)
pd = fitdist(spec19(61,rec19 <=190)', 'gamma');
g1 = pdf(pd,  0:0.1:5);
hold on
plot( 0:0.1:5, g1, 'b','LineWidth', 2)
[h2, x2 ] =  hist(spec19(61,rec19 >190), x1);

bh2 = bar(x2, h2/sum(h2)/(x2(2)-x2(1)), 1)
set(bh2,'facecolor','g');
set(get(bh2,'Children'),'FaceAlpha',0.2)
pd = fitdist(spec19(61,rec19 >190)', 'gamma');
g2 = pdf(pd, 0:0.1:5);
plot( 0:0.1:5, g2, 'g','LineWidth', 2)
xlabel('Magnitude', 'fontsize',fs); ylabel('Histogram frequency', 'fontsize',fs)
hold off


subplot(2,1,1)
plot((10:128)*48/128, spec6(10:end, 20:22))
xlim([4, 48]); xlabel('Frequency (Hz)', 'fontsize',fs); ylabel('Magnitude', 'fontsize',fs);
title('a', 'fontsize',fs+2)
subplot(2,1,2)
a = spec6(72,:);
a = a(a~=0)';
pd = fitdist(a, 'Gamma');
[h1, x1] = hist(a, 20);
g1 = pdf(pd, x1);
bar(x1, h1/sum(h1)/(x1(2)-x1(1)), 1); %xlim([0,4])
hold on
plot(x1, g1, 'r')
hold off
xlabel('Magnitude', 'fontsize',fs); ylabel('Histogram frequency', 'fontsize',fs);
title('b', 'fontsize',fs+2)

imagesc(log(spec6(10:end, :))); caxis(prctile(log(spec6(:)), [10, 95])); set(gca,'YDir','normal'); ylabel('frequency')

a = spec6(50,:);
a = a(a~=0)';
pd = fitdist(a, 'Gamma');
[h1, x1] = hist(a, 20);
g1 = pdf(pd, x1);
%subtightplot(8,8,fi, [0.02, 0.01])

bar(x1, h1/sum(h1)/(x1(2)-x1(1)), 1); %xlim([0,4])
hold on
plot(x1, g1, 'r')
hold off


%play
for fi = 1:64
    f = 2*fi;
    a = spec19(f, rec19 <=190);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h1, x1] = hist(a, 20);
      g1 = pdf(pd, x1);
    subtightplot(8,8,fi, [0.02, 0.01])
     
    bar(x1, h1/sum(h1)/(x1(2)-x1(1)), 1); %xlim([0,4])
    hold on
    plot(x1, g1, 'r')
    hold off
      
        a = spec19(f, rec19 >190);
    a = a(a~=0)';
    pd = fitdist(a, 'Gamma');
      [h2, x2] = hist(a, 10);
      g2 = pdf(pd, x2);
    %subtightplot(8,8,fi, [0.02, 0.01])
     bar(x2, h2/sum(h2)/(x2(2)-x2(1)), 1); %xlim([0,4])
    hold on
    plot(x2, g2, 'r')
    hold off
end
