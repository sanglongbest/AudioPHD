%GM_PHD_Simulate_Plot
%Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au 
%error_ellipse is by AJ Johnson, taken from Matlab Central http://www.mathworks.com.au/matlabcentral/fileexchange/4705-errorellipse

%This file plots the current measurements, true target position and estimated target
%position, as well as a history of estimated target positions.

disp('Plotting...')
%Measurements
figure(1);
if(PLOT_ALL_MEASUREMENTS == false)
    clf;%Only plot the most recent measurement.
end
hold on;

%Plot all measurements, including clutter, as black 'x'
if(~isempty(Z))
    plot(k,Z(1,:), 'xk');
end
%Plot noisy measurements of true target position(s), as black 'o'
% if(~isempty(zTrue))
%     plot(zTrue(1,:), zTrue(2,:), 'ok');
% end

%Plot target 1 true position as red dots
plot(0:k,simTarget1History(1,:), '.-r');
%Plot target 2 true position as blue dots
plot(0:k,simTarget2History(1,:), '.-b');

%Plot tracked targets as magenta circles
if(~isempty(X_k_history))
    plot(k,X_k_history(1,:), 'om');
end
xlabel('time');
ylabel('tau');
title('Simulated targets and measurements');
axis square;

%For extracted targets, plot latest target(s) as cyan triangle, and draw an
%error ellipse to show uncertainty
if(~isempty(X_k))
    plot( k, X_k(1,:), '^c');
    [nRows, nCols] = size(X_k);
    for c = 1:nCols
       thisMu = X_k(1:2, c);
       covRange = calculateDataRange4(c);
       thisCov = X_k_P(:,covRange);
       thisCov = thisCov(1:2, 1:2); %We only care about position
       error_ellipse(thisCov, thisMu);
    end
    if(DRAW_VELOCITIES == 1)
      %Draw velocities of targets   
      quiver(X_k(1,:), X_k(2,:), X_k(3,:), X_k(4,:))          
    end
end

%Individual X and Y components of measurements
figure(2);
subplot(2,1,1);
plot(k, Z(1,:), 'xk');%X coord of clutter measurements
if(~isempty(simTarget1History))
     plot(k, simTarget1History(1,:), 'or');
end
plot(k, simTarget1History(1,:), 'xk');%X coord of true measurement

if(~isempty(X_k))
    plot(k, X_k(1,:), 'om');
end
subplot(2,1,2);
plot(k, Z(1,:), 'xk');%Y coord of clutter measurements
if(~isempty(simTarget2History))
     plot(k, simTarget2History(1,:), 'ob');
end
if(~isempty(X_k))
    plot(k, X_k(2,:), 'om');
end

%OSPA error metric plot
if(CALCULATE_OSPA_METRIC == 1)
    figure(3);
    clf;
    hold on;
    plot(metric_history, 'x-b');

    axis([0 100 0 cutoff_c]);
    xlabel('Simulation step');
    ylabel('OSPA error metric (higher is worse)');
    title('OSPA error metric for this test');
end