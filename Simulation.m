
function[EZ1,EZ2,Epower1,Epower2,runningEZ1,runningEZ2] = Simulation(delta,lambda,mu);
% Simulation for case study product form queue

% n = 3;
% v0 = 1;
% edges = [1 2; 2 3];
% r = [0.01 0.005];
% x = [0 0];
% Ymat = zeros(n,n);
% for i=1:length(edges)
%     Ymat(edges(i,1),edges(i,2)) = -1/(r(i)+1i*x(i));
%     Ymat(edges(i,2),edges(i,1)) = -1/(r(i)+1i*x(i));
% end
% for i=1:n
%     Ymat(i,i) = -sum(Ymat(i,:));
% end
% Zmat = inv(Ymat(2:end,2:end));
% Rmat = real(Zmat);
% Xmat = imag(Zmat);
% Amat= [1/((1-delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat   -1/((1+delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat;
%     1/((1-delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat   -1/((1+delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat;
%     -1/((1+delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat    1/((1-delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat;
%     -1/((1+delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat    1/((1-delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat];
% r1 = .01;
% r2 = .015;


powerallocation = load(['LambdaMatrices/LambdaMatrixSimulation' num2str(delta) '.mat']);


arrivals = [exprnd(1/(2*lambda))];
chargesnode1 = [];
chargesnode2 = [];

runtime = 10000;
warmuptime = 0.95*runtime;
t = 0;
previous_t = 0;
power1 = 1;
power2 = 1;
sum_EVs_node1 = 0;
sum_EVs_node2 = 0;
sum_power_node1 = 0;
sum_power_node2 = 0;
runningEZ1 = [0];
runningEZ2 = [0];
while t < runtime
    previous_t = t;
    if isempty(chargesnode1)
        if isempty(chargesnode2)
            [t, event_type] = min([arrivals(1)]);
        else
            [t, event_type] = min([arrivals(1),10^5,chargesnode2(1)]);
        end
    else
        if isempty(chargesnode2)
            [t, event_type] = min([arrivals(1),chargesnode1(1)]);
        else
            [t, event_type] = min([arrivals(1),chargesnode1(1),chargesnode2(1)]);
        end
    end
%     disp("delta = " + num2str(delta) + " t = "+num2str(t) + " type = " + num2str(event_type))
    % save performance measures
    if t > warmuptime % exclude warm-up time
        sum_EVs_node1 = sum_EVs_node1 + (t-previous_t)*length(chargesnode1);
        sum_EVs_node2 = sum_EVs_node2 + (t-previous_t)*length(chargesnode2);
        sum_power_node1 = sum_power_node1 + (t-previous_t)*power1;
        sum_power_node2 = sum_power_node2 + (t-previous_t)*power2;
        runningEZ1(end+1) = sum_EVs_node1/(t-warmuptime);
        runningEZ2(end+1) = sum_EVs_node2/(t-warmuptime);
%         runningEpower1 = sum_power_node1/(t-warmuptime);
%         runningEpower2 = sum_power_node2/(t-warmuptime);
    end
    
    % handle event
    if event_type == 1 % arrival 
        arrivalnode = binornd(1,0.5)+1;
        if arrivalnode == 1
            chargesnode1(end+1) = t+exprnd(mu)/power1;
%             disp("dep planned node1 at " + chargesnode1(end) + " power1 " + power1)
        else
            chargesnode2(end+1) = t+exprnd(mu)/power2;
%             disp("dep planned node2 at " + chargesnode2(end) + " power2 " + power2)
        end
        arrivals(end+1) = t+exprnd(1/(2*lambda)); % plan new arrival
        arrivals = arrivals(2:end); % delete event from eventq
        
    elseif event_type == 2 % charged node 1
        if length(chargesnode1) > 1
            chargesnode1 = chargesnode1(2:end);
        else
            chargesnode1 = [];
        end
    elseif event_type == 3 % charged node 2
        if length(chargesnode2) > 1
            chargesnode2 = chargesnode2(2:end);
        else
            chargesnode2 = [];
        end
    end
    
    % update power allocation
    previous_power1 = power1;
    previous_power2 = power2;
    z = [length(chargesnode1) length(chargesnode2)];
%     [alc_mod,alc_socp,alc_ldf] = OPF_SOCP_and_R_and_LDF(n-1,Ymat,Zmat,Amat,v0,delta,[1 1],0,0,z);
    if z(1) > 0
        power1 = powerallocation.socp(z(1)+1,z(2)+1,2);
%         power1 = alc_socp(2);
% power1 = ((1-delta)*delta)/(r1*(z(1)+z(2)));%(1-(1-delta)^2)/(2*r1*(z(1)+z(2)));
    else
        power1 = 1;
    end
    if z(2) > 0
        power2 = powerallocation.socp(z(1)+1,z(2)+1,3);
%         power2 = alc_socp(3);
% power2 = ((1-delta)*delta)/(r2*(z(1)+z(2))); %(1-(1-delta)^2)/(2*r2*(z(1)+z(2)));
    else
        power2 = 1;
    end
    chargesnode1 = t*ones(1,z(1)) + (chargesnode1 - t*ones(1,z(1))).*(previous_power1/power1);
    chargesnode2 = t*ones(1,z(2)) + (chargesnode2 - t*ones(1,z(2))).*(previous_power2/power2);     
end

EZ1 = sum_EVs_node1/(t - warmuptime);
EZ2 = sum_EVs_node2/(t - warmuptime);
Epower1 = sum_power_node1/(t - warmuptime);
Epower2 = sum_power_node2/(t - warmuptime);
end