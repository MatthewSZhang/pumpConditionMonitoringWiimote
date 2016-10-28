function X_hp = filterHpfFirUsingLpfMaPhaseCorrect(X,plotOption)
% HPF = Signal - LPF(moving average)

% Moving average filter coeff
% Ad-hoc: Currently, fixed at 4th order!
a = 1;
b = [1/4 1/4 1/4 1/4];
% AM Modified 25 Oct 2016: Add delay to take care of phase shift
% Take care of phase delay by appending zeros to original signal
% X_ma = filter(b,a,X);
D = 2; % D = filter order/2
X_ma = filter(b,a,cat(1,X,zeros(D,1)));
X_ma = X_ma(D+1:end);
X_hp = X-X_ma;        

if (plotOption)        
    timeStamp = (1:length(X))';
    yMax = max(cat(1,X,X_ma,X_hp));
    yMin = min(cat(1,X,X_ma,X_hp));  

    % Plot signal and low-passed
    figure(1);
    subplot(2,1,1);
    plot(timeStamp,X,'b');
    hold on;            
    plot(timeStamp,X_ma,'r');
    ylabel('Y Low-passed');
    axis([0 length(X) yMin yMax]);
    hold off;
    legend('X','Xma');
    grid on; 

    figure(1);
    subplot(2,1,2);
    plot(timeStamp,X,'b');
    hold on;  
    plot(timeStamp,X_hp,'r');
    hold off;
    axis([0 length(X) yMin yMax]);
    ylabel('Y High-passed');
    legend('X','Xhp');  
    grid on; 

%             % Plot diff between filtfilt and HPF=LPF-MA            
%             % Moving average filter coeff   
%             figure(fileId+100);subplot(3,1,2); 
%             plot(timeStamp,X_hp,'b');
%             hold on;            
%             plot(timeStamp,X_hpma,'r');
%             axis([0 yMax -.4 .4]);
%             hold off;
%             legend('Yhp','Yhpma');
% 
%             figure(fileId+100);subplot(3,1,3);
%             plot(timeStamp,X_hp-X_hpma);
%             axis([0 yMax -.4 .4]);            
%             legend('Yhp-Yhpma');
%             suptitle('Yhp vs. Yhpma');
    drawnow;
end