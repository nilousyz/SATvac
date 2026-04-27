function [M] = SplitTransitionM(M,t,steptime,stoptime,pdivM,pleaveM)

if ~isempty(M.Transit)
    idx = M.Transit(1,:) < t+0.1;
    M.CellCount = M.CellCount - sum(M.Transit(2,idx));          
    M.ToAdd = M.ToAdd + (2*sum(M.Transit(2,idx)));              
    if M.ToAdd > 5300000
        warning('excessive M')
        error('Error Occurred')
    end
    M.Transit = M.Transit(:,~idx);                              
end


if  M.ToAdd > 0
    M.CellCount = M.CellCount + M.ToAdd;
    AgeTransit = (steptime.*geornd(pdivM,1,M.ToAdd)) + t + steptime;
    AgeLeaving =  (steptime.*geornd(pleaveM,1,M.ToAdd)) + t + steptime;
    idx = AgeTransit >= AgeLeaving;
    AgeTransit(idx | AgeTransit > stoptime) = [];
    AgeLeaving(~idx | AgeLeaving > stoptime) = [];
    M.Transit = UpdateCellVectors(M.Transit, AgeTransit,1);  
    M.Leaving = UpdateCellVectors(M.Leaving, AgeLeaving,1);  
    M.ToAdd = 0 ;
end


end