function [Ef,rate] = SplitTransitionEf(Ef,t,steptime,stoptime,pdiv,pleave,rate)

if ~isempty(Ef.Transit)
    idx = Ef.Transit(1,:) < t+0.1;
    Ef.CellCount = Ef.CellCount - sum(Ef.Transit(2,idx));
    Ef.ToAdd = Ef.ToAdd + (2*sum(Ef.Transit(2,idx)));
    rate = (Ef.ToAdd/Ef.CellCount)*100;
    if Ef.ToAdd > 1000000000
        warning('Excessive Ef')
        error('Error Occurred')
    end
    Ef.Transit = Ef.Transit(:,~idx);

end

if Ef.ToAdd > 0
    Ef.MToAdd = round((rand(1,1)*0.06)* Ef.ToAdd);
    Ef.CellCount = Ef.CellCount + Ef.ToAdd - Ef.MToAdd;
    Ef.ToAdd = Ef.ToAdd - Ef.MToAdd;
    AgeTransit = (steptime.*geornd(pdiv,1,Ef.ToAdd)) + t + steptime;
    AgeLeaving =  (steptime.*geornd(pleave,1,Ef.ToAdd)) + t + steptime;
    idx = AgeTransit >= AgeLeaving;
    AgeTransit(idx | AgeTransit > stoptime) = [];
    AgeLeaving(~idx | AgeLeaving > stoptime) = [];
    Ef.Transit = UpdateCellVectors(Ef.Transit, AgeTransit,1);
    Ef.Leaving = UpdateCellVectors(Ef.Leaving, AgeLeaving,1);
    Ef.MToAdd = 0;
    Ef.ToAdd = 0 ;
end

end