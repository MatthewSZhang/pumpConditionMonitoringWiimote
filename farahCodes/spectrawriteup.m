clear
close all
% Condition monitoring

%General gist: 

%Design of high pass filte:Done fairly arbitrarily, based on what looks OK
Fs = 96;  % Sampling Frequency

Fstop = 2;               % Stopband Frequency
Fpass = 4;               % Passband Frequency
Dstop = 0.0001;          % Stopband Attenuation
Dpass = 0.057501127785;  % Passband Ripple
dens  = 20;              % Density Factor
% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop, Fpass]/(Fs/2), [0 1], [Dstop, Dpass]);
% Calculate the coefficients using the FIRPM function.
Highpass  = firpm(N, Fo, Ao, W, {dens});

LFpass = 3; 
LFstop = 5; 
LDpass = 0.057501127785;  
LDstop = 0.0001; 
Ldens  = 20;   
[N, Fo, Ao, W] = firpmord([LFpass, LFstop]/(Fs/2), [1 0], [LDpass, LDstop]); 
Lowpass  = firpm(N, Fo, Ao, W, {Ldens});

% AM Modified
figure(1);subplot(2,1,1);plot(Lowpass);title('LowPass');
figure(1);subplot(2,1,2);plot(Highpass);title('HighPass');

%Size of windows and their overlap. Trade off between frequency/time
%specifity and smoothness
cliplength = 256; overlap = 128;

% AM Modified
rawDataPath = 'C:\Users\engs1602\research\data\Data Feb-Mar 2016\Raw Data - wiimote';

% Begin code
spectraThisPump = [];
idx = 1;
% AM Modified
for rec = 1:16%1:313
    %open each recording - should really be in seperate function!
    % AM Modified
%     file = fopen(strcat(num2str(rec), '.txt'), 'rt');
    file = fopen(fullfile(rawDataPath,sprintf('%d.txt',rec)),'rt');
    raw = textscan(file, '%f:%f:%f:%f \t %f \t %f \t %f \t %f \t %f \t');
    fclose(file);
    time = raw{4}/1000 + raw{3} + 60*raw{2};
    Y = raw{6}; 
    Z = raw{7};
    Y = Y - mean(Y);
    Z = Z - mean(Z);
    
    if length(Y) > 396
        % lowpass filter and arc calculation for gross movement
        Y_s = filtfilt(Lowpass,1,Y);
        Z_s = filtfilt(Lowpass,1,Z);
        arc = atand(Y_s./(Z_s+Z));

        % AM Modified
        figure(rec);subplot(4,1,1);plot(time,Y);ylabel('Y');grid on;
        figure(rec);subplot(4,1,2);plot(time,Y_s);ylabel('Y Low-passed');grid on;
        
        % Just removes stationary periods at start and end of recording
        % How to pick the threshold? Very ad-hoc!
        threshold = 0.25;
        interval2 = sortoutpeaks3FarahNA(Y_s, threshold); 
        time = time(interval2);
        Y = Y(interval2);
        % AM Modified
        figure(rec);subplot(4,1,3);plot(time,Y);ylabel('Y ends removed');grid on;  
        
         if length(Y) > 396
            intstart = 1;
            %Yf = filter(Highpass, Y);
            % Apply highpass filter to whole recording
            Yf = filtfilt(Highpass, 1, Y);
            % AM Modified            
            figure(rec);subplot(4,1,3);plot(time,Yf);ylabel('Y High-passed');

            %for each window in the clip
            while intstart < length(interval2) - cliplength
                intend = intstart + cliplength;
                %apply fft to the window
                Yfhat = fft(Yf(intstart:intend), cliplength);
                % save this up to the nyquist limit. 
                spectra(:,idx) = abs(Yfhat(1:cliplength/2));
                rec_no(idx) = rec;
                % calculate amount of movemnet in the window
                spd(idx) = sum(abs(diff(arc(intstart:intend))));
                idx = idx+1;
                intstart = intend;
            end
            figure(100+rec);
            subplot(3,1,1:2)
            imagesc(spectra(10:end, :)); caxis([0, 2.5])
%             set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'Fontsize', 12);
%             set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'Fontsize', 12);
            set(gca,'YDir','normal'); ylabel('frequency')
            subplot(3,1,3)
            plot(spd); xlim([1,size(spectra,2)]); %ylim([150, 230])
%             set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'Fontsize', 12);
%             set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'Fontsize', 12);
            xlabel('time'); ylabel('speed')
         end
    end    
end

% % AM Modified
% %looking at it
% imagesc(lspectra(10:end, loc_num(rec_no) == l))

% % AM Modified
% figure(6);
% subplot(3,1,1:2)
% plot(spectra(:, 10:15), 'b')
% hold on
% plot(median(spectra(:, :), 2), 'LineWidth', 3)
% hold off
% plot(spectra(:, 10:15), 'b')   
% xlabel('Frequency'); ylabel('Magnitude')
   
% figure(7)
% subplot(3,1,1:2)
% imagesc(spectra(10:end, :)); caxis([0, 2.5])
% set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'Fontsize', 12);
% set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'Fontsize', 12);
% set(gca,'YDir','normal'); ylabel('frequency')
% subplot(3,1,3)
% plot(spd); %xlim([1,73]); ylim([150, 230])
% set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'Fontsize', 12);
% set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'Fontsize', 12);
% xlabel('time'); ylabel('speed')

% l= 1;
% for f = 1:64
%     a = spectra(f, loc_num16(rec_no) == l);
%     a = a(a~=0)';
%     pd = fitdist(a, 'Gamma');
%       [h1, x1] = hist(a, 30);
%       g1 = pdf(pd, x1);
%     subtightplot(8,8,f, [0.02, 0.01])
%       bar(x1, g1, 'g')      
% end