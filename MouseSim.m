
function [dat,total] = MouseSim(NA0, NP0, pAPbinding, IntroductionRate, BPLeavingRate, pEHleaving, stoptime, steptime, X7 , X8 , X9 , X10 , X11)

pJ2J.Mean = X7;
pJ2J.STD = X8;
pJ2J.Scaling = 1;
pIJleaving = X9;
pSplit.Mean = X10;
pSplit.STD = X11;
pSplit.Scaling = 1;

% Set up the variables we care about
dat = zeros(1, 19);
dat(1, 2) = 0;
dat(1,12) = 0;
total = zeros(1,2);
i = 1; % starting at row 1
bind_count = [];
rate= 0;
Ef.ToAddFromJ = 0;

% List of cell types
cellTypes = {'A','K','C','D','E','F','G','H','I','J','P','AP','BP','DP','M','Ef'};

for i = 1:numel(cellTypes)
    eval([cellTypes{i} '.CellCount = 0;']);
    eval([cellTypes{i} '.ToAdd = 0;']);
    eval([cellTypes{i} '.Transit = [];']);
    eval([cellTypes{i} '.Leaving = [];']);
end

% vector of variables in various states, each holds the age of each cell

P.AgeBinding = [];
P.AgeBindingM = [];
P.AgeBindingEf = [];
P.AgeBindingA = [];

pAP.Mean = 0.458;
pAP.STD = 0.15;
pAP.Scaling = 0.25;
pAP.Agemax = 1.1;
pAP.Agemin = 0;
pAP.Shift = 0;

pBP.Mean = 21;
pBP.STD = 3;
pBP.Scaling = 1;

TotalAadded = 0;
TotalPadded = 0;
TotalPremoved = 0;
K.AgeMigrate = [];
Blood.CellCount = 0;

% for each timestep of the simulation
t = 0;
while t <= stoptime
    %% Division and leaving rates for effector and memory

    kdiv = 0.015;
    adiv = 0.15;

    kleave= 0.0078;
    aleave = 0.03;

    pdiv = adiv * exp(-kdiv * t);
    pleave = aleave * exp(-kleave * t);

    kdivM = 0.014;
    adivM = 0.1;

    kleaveM= 0.0078;
    aleaveM = 0.025;

    pdivM = adivM * exp(-kdivM * t);
    pleaveM = aleaveM * exp(-kleaveM * t);

    %% A and P Introduction

    if (t < IntroductionRate)
        if (TotalAadded < NA0)
            [A.ToAdd] = CellIntroduction (2,IntroductionRate,t,steptime,NA0,TotalAadded);
            TotalAadded = TotalAadded + A.ToAdd;
        end
    end
    secondInjTime=24;
    IntroductionRateP=24;
    if (t < (secondInjTime+IntroductionRateP))
        if (TotalPadded < NP0 && t > secondInjTime)
            [P.ToAdd] = CellIntroduction (2,IntroductionRateP,t-secondInjTime,steptime,NP0,TotalPadded);
            TotalPadded = TotalPadded + P.ToAdd;
        end
    end

    if (t > IntroductionRate && BP.CellCount > 0 && TotalPremoved < NP0 && BPLeavingRate > 0)
        TotalPThisTimestep = round(BPLeavingRate * (t - IntroductionRate + steptime) * NP0 / IntroductionRate);
        P.ToRemove = min(max([TotalPThisTimestep - TotalPremoved, 0]), BP.CellCount); 

        if (P.ToRemove > 0)
            TotalPremoved = TotalPremoved + P.ToRemove;
            if (P.ToRemove >= BP.CellCount)
                BP.CellCount = 0;
                BP.Transit=[];
            else 
                for PIndex = 1:P.ToRemove
                    idx = randi([1 size(BP.Transit,2)], 1);
                    BP.Transit = UpdateCellVectors(BP.Transit, BP.Transit(1,idx),-1); 
                    BP.CellCount = BP.CellCount-1;
                end
            end
        end
    end

    %% Cell Leaving to K

    % for Ef
    [Ef,K] = CellLeavingToNewK (Ef,K,t,steptime);
    % for M
    [M,K] = CellLeavingToNewK (M,K,t,steptime);
    % for I
    [I,K] = CellLeavingToNewK (I,K,t,steptime);
    % for J
    [J,K] = CellLeavingToNewK (J,K,t,steptime);
    % for E
    [E,K] = CellLeavingToNewK (E,K,t,steptime);
    % for F
    [F,K] = CellLeavingToNewK (F,K,t,steptime);
    % for G
    [G,K] = CellLeavingToNewK (G,K,t,steptime);
    % for H
    [H,K] = CellLeavingToNewK (H,K,t,steptime);

    %% Cell Transition

    % M --> 2M transition
    [M] = SplitTransitionM(M,t,steptime,stoptime,pdivM,pleaveM);
    % Ef --> 2Ef transition
    [Ef,rate] = SplitTransitionEf(Ef,t,steptime,stoptime,pdiv,pleave,rate);
    % J --> invisible  or J --> 2J
    [J,M,Ef] = SplitTransitionJ(J,pJ2J,pIJleaving,t,steptime,stoptime,M,Ef,pdiv,pleave,pdivM,pleaveM);
    % I -> 2J transition
    [I,J] = SplitTransitionI(I,J,t);
    % H -> 2I transition
    [H,I] = SplitTransition(H,I,pSplit,pEHleaving,t,steptime,stoptime);
    % G -> 2H transition
    [G,H] = SplitTransition(G,H,pSplit,pEHleaving,t,steptime,stoptime);
    % F -> 2G transition
    [F,G] = SplitTransition(F,G,pSplit,pEHleaving,t,steptime,stoptime);
    % E -> 2F transition
    [E,F] = SplitTransition(E,F,pSplit,pEHleaving,t,steptime,stoptime);
    % D -> 2E transition
    [D,E] = SplitTransitionD(D,E,pSplit,pEHleaving,t,steptime,stoptime);
    % C -> 2D transition
    [C,D] = SplitTransitionC(C,D,pSplit,t,steptime,stoptime);

    %% BP Transition

    [BP,C,P] = BPTransition(BP,C,pSplit,t,P,steptime,stoptime);

    %% AP Transition

    [AP,BP] = APTransition(AP,BP,pBP,t,steptime);

    %% A + P binding   A + P --> AP transition using normal distributions

    if (A.ToAdd > 0)
        A.CellCount = A.CellCount + A.ToAdd;
        A.ToAdd=0;
    end

    if (P.ToAdd > 0)
        P.CellCount = P.CellCount + P.ToAdd;
        Prob = pAPbinding*(A.CellCount / (A.CellCount + P.CellCount));
        P.AgeBinding = [P.AgeBinding (steptime.*geoinv(rand(1,P.ToAdd),Prob) + t )];
        P.ToAdd=0;
    end

    if (A.CellCount > 0 && P.CellCount > 0)
        idx = find(P.AgeBinding <= t);
        if size(idx,2)> A.CellCount
            idx=idx(1:A.CellCount);
        end
        A.CellCount = A.CellCount - size(idx,2);
        P.CellCount = P.CellCount - size(idx,2);
        bind_count(i) = size(idx,2);
        P.AgeBinding(idx) = [];
        AP.ToAdd = size(idx,2);

    end

    if AP.ToAdd > 0
        AgeTransit = (inverseCDF(AP.ToAdd,pAP) + t + steptime);
        idx = AgeTransit > stoptime;
        AgeTransit(idx) = [];
        AP.Transit = UpdateCellVectors(AP.Transit, AgeTransit,1);  
        AP.CellCount = AP.CellCount + AP.ToAdd;
        AP.ToAdd=0;
    end
    i = i + 1; % next iteration
    t = round(t + steptime,1); % next time step

    %% store all of the relevant data
    dat(i,:) = [t,A.CellCount,K.CellCount,C.CellCount,D.CellCount,E.CellCount,F.CellCount,G.CellCount,H.CellCount,I.CellCount,J.CellCount,P.CellCount,AP.CellCount,BP.CellCount,DP.CellCount,M.CellCount,Ef.CellCount,Blood.CellCount,rate];
    total(i,:) = [t,(sum(dat(i,[2,4:11,13:17])))];
end

end