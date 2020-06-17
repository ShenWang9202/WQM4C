 function [UeachMinforEPANET,U_C_B_eachStep,PreviousSystemDynamicMatrix] = ObtainControlAction(CurrentValue,IndexInVar,aux,ElementCount,q_B,x_estimated,PreviousValue)
% Here we need to know the current concentration vector x which includes
% the c^J, c^R, c^TK, c^P, c^M, and c^W.

delta_t = CurrentValue.delta_t;

JunctionCount= ElementCount.JunctionCount;
ReservoirCount = ElementCount.ReservoirCount;
TankCount = ElementCount.TankCount;

JunctionCount = double(JunctionCount);
ReservoirCount = double(ReservoirCount);
TankCount = double(TankCount);
nodeCount = JunctionCount + ReservoirCount + TankCount;

flowRate_B = aux.flowRate_B;
q_B = aux.q_B;
Price_B = aux.Price_B;

Np = CurrentValue.Np;
[A,B,C] = ObtainDynamicNew(CurrentValue,IndexInVar,aux,ElementCount,q_B);

PreviousSystemDynamicMatrix = struct('A',A,...
    'B',B,...
    'C',C);


% MPC data
[W,Z,Ca,PhiA,GammaA] = MPCScale(A,B,C,Np);
% size(W)
% size(Z)
% size(Ca)
% size(PhiA)
% size(GammaA)
% full(W)
% full(Z)
% full(Ca)
% full(PhiA)
% full(GammaA)

% Xa = [Delta x;y]
x = x_estimated; 
[nx,~] = size(x_estimated);
yk = C*x;
[ny,~] = size(yk);
Deltax = zeros(nx,1);
Xa = [Deltax;yk];
% only select the booster location;
nu = nodeCount;
n_deltau = Np*nu;
%R = 1000*speye(n_deltau,n_deltau); Change this R to change the smoothness
%of control action.
R =  speye(n_deltau,n_deltau);
%R =  sparse(n_deltau,n_deltau);
% for i = 1:Np
%     R_i = eye(nu,nu);
% %     R_i = zeros(nu,nu);
% %     for j = 1:JunctionCount
% %         R_i(j,j) = 1; 
% %     end
%     R = blkdiag(R, R_i);
% end
% reference
n_ref = ny*Np;
reference = Constants4Concentration.reference;
reference = reference*ones(n_ref,1);
% only select the pipes now.
n_Q = Np*ny;
Q = speye(n_Q,n_Q);
% for i = 1:Np
%     Q_i = eye(ny,ny);
% %    Q_i = zeros(ny,ny);
% %     for j = (nodeCount + 1):(nodeCount + NumberofSegment*PipeCount)
% %         Q_i(j,j) = 1; 
% %     end
%     Q = blkdiag(Q, Q_i);
% end
% This is just C_B;

%DeltaU1 = (R + Z'*Q*Z)^(-1)*Z'*Q*(reference-W*Xa);

b = [];
%b = ones(n_deltau,1);
for i = 1:Np
    b = [b; (Np-i)*Price_B];
end
Price_Weight = Constants4Concentration.Price_Weight;

%DeltaU = (R + Z'*Q*Z)^(-1)*(Z'*Q*(reference-W*Xa)-Price_Weight*b);
DeltaU = (R + Z'*Q*Z)\(Z'*Q*(reference-W*Xa)-Price_Weight*b);
% We reshape it based on the number of nodes.
[m_Delta,~] = size(DeltaU);
column = m_Delta/nu;
DeltaU = reshape(DeltaU,nu,column);
% % Just pick the Junction2 out
% DeltaU_Junction2 = DeltaU(1,:);
% Times the flowRate_B to get the mass which is in mg/min
%DeltaU =  DeltaU .* q_B .* Constants4Concentration.Gallon2Liter;
% Accumulate Delta U and make it as U;

% Note DeltaU is a relative value in terms of each step (delta_t)
[m_Delta,n_Delta] = size(DeltaU);
U = zeros(m_Delta,n_Delta);
U(:,1) = DeltaU(:,1);
for i = 1:n_Delta-1
    U(:,i+1) = U(:,i) + DeltaU(:,i+1);
end

% Note that this U is still a relative value (still a delta value) in terms of its previous value
% (5 minutes ago)

% Now the U_C_B is in mg/L.
U_C_B = U;

U_C_B_eachStep = U_C_B;
% Next, we obtain the true U
% Get the prevoius U;
PreviousU_C_B_eachStep = PreviousValue.U_C_B_eachStep;
% Get the Current asolute U;
U_C_B_eachStep = U_C_B_eachStep + PreviousU_C_B_eachStep(:,end);


% make sure the U_C_B_eachStep IS NOT NEGATIVE
% Because we only can add concentration to it, and cannot extract
% concnetration from water.
[m_UCB,n_UCB] = size(U_C_B_eachStep);
for i = 1:m_UCB
    for j = 1:n_UCB
        if(U_C_B_eachStep(i,j) < 0)
            disp('U_C_B_eachStep negative')
            U_C_B_eachStep(i,j) = 0;
        end
    end
end


% fix this, since we need to calcuate the average value for EPAENT in a
% minute. Our LDE model can apply in second, but due to the simulation step
% for EPANET is one minute, we only can apply control action at the integer
% minuter (1 2 3 4 .... minutes), but for LDE it is (1.2 1.3 1.4 minutes)

% even though the LDE and EPANET are close superly, the control action
% applied to the epanet and lde is different. This makes the final control
% effect is a little bit different.


% Next, we calcuate the realtive U for EPANET.



% This U is acutally C_B; UforEPANET is the source quality for EPAENT
% software since it requirs mg/min

% What we get here U (mg/L) is the concentration for each step,and our each step is in second instead of minute. 
% Besides that, we assume the booster flow rate is q_B (gallon/min),so in order to get the mg/sec, we should convert
% flow reate to L/sec first,

% Now the UforEPANET is in mg for each step).
UforEPANET =  U_C_B_eachStep .* q_B .* Constants4Concentration.Gallon2Liter ./ Constants4Concentration.MinInSecond * delta_t;


% Get the prevoius U;
PreviousUforEPANET = PreviousValue.UeachMinforEPANET;

% Get the absolute U first
%UforEPANET = UforEPANET + PreviousUforEPANET(:,end);

% now this is the U for next Hq_min = 5 mins here
% But our Quality Time Step is set to 1 min, so we need to pick out the
% value at each minute.
IndexofUeachMin = [0];
Hq_min = Constants4Concentration.Hq_min;

for i = 1:Hq_min
    IndexofUeachMin = [IndexofUeachMin round(i * Constants4Concentration.MinInSecond/delta_t)];
end

UeachMinforEPANET = [];
for i = 1:Hq_min
    startIndex = IndexofUeachMin(i) + 1;
    endIndex =  IndexofUeachMin(i+1);
    %startIndex:endIndex
    eachMinMass = sum(UforEPANET(:,startIndex:endIndex),2);
    UeachMinforEPANET = [UeachMinforEPANET eachMinMass];
end

%UeachMinforEPANET = UeachMinforEPANET + PreviousUforEPANET(:,end);

% make sure the UeachMinforEPANET IS NOT NEGATIVE
[m_UCB,n_UCB] = size(UeachMinforEPANET);
for i = 1:m_UCB
    for j = 1:n_UCB
        if(UeachMinforEPANET(i,j) < 0)
            disp('U_C_B_eachMin negative')
            UeachMinforEPANET(i,j) = 0;
        end
    end
end

%U_C_B_eachStep
%UeachMinforEPANET

