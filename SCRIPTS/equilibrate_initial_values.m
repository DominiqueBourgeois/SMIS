function init=equilibrate_initial_values(K)
% Get steady state values for a kinetic system 
% K(S1,S2) is the exchange rate between state S1 and S2
% Returns the fraction of molecules in each state at t=0;

%Modified DB February 20200: adjust break_value

if isempty(K); init=[]; return; end

s=size(K,1);
init=1/s*ones(1,s); % Initialize init
% break_val=1e-6; % Stop when change in C become smaller than that
precision_level=1e-3;

max_k=max(K(:));
if max_k==0; return; end

min_k=min(min(K(K>0)));
dt=0.1/max_k; 
T=10/min_k; 
N=fix(T/dt); % # of steps

C=init; % Initial concentrations
dC=zeros(1,s); 
for i=1:N
    for j=1:s
       dC(j)=-sum(C(j)*K(j,:))*dt+ sum(C'.*K(:,j))*dt ;
    end
    C=C+dC;
    
    if i==1
      break_val=max([max(abs(dC))*precision_level,1e-10]); % Break value should adapt to chosen dt and thus dC
    end
    
    if max(abs(dC))<break_val; break; end
end
init=C;
end
