function [Cell,M,Ef] = SplitTransitionJ(Cell,pSplit,CellLeavingProb,t,steptime,stoptime,M,Ef,pdiv,pleave,pdivM,pleaveM)

if ~isempty(Cell.Transit)
    idx = Cell.Transit(1,:) < t+0.1;
    Cell.CellCount = Cell.CellCount - sum(Cell.Transit(2,idx));         
    Cell.ToAdd = Cell.ToAdd + (2*sum(Cell.Transit(2,idx)));             
    Cell.Transit = Cell.Transit(:,~idx);                                
end

if  Cell.ToAdd > 0
    Ef.ToAddFromJ = round (binornd(Cell.ToAdd, 0.85));
    Cell.ToAdd = Cell.ToAdd - Ef.ToAddFromJ ;
    M.ToAdd = round(binornd(Cell.ToAdd, 0.40));
    Cell.ToAdd = Cell.ToAdd - M.ToAdd ;
    Ef.ToAdd = Ef.ToAdd + Ef.ToAddFromJ;
    AgeLeaving = (steptime.*geoinv(rand(1,Cell.ToAdd),CellLeavingProb) + t + steptime );
    AgeTransit = (inverseCDF_j(Cell.ToAdd,pSplit) + t + steptime );
    idx = AgeTransit >= AgeLeaving;
    AgeTransit(idx | AgeTransit>stoptime) = [];
    AgeLeaving(~idx | AgeLeaving>stoptime) = [];
    Cell.Transit = UpdateCellVectors(Cell.Transit, AgeTransit,1);  
    Cell.Leaving = UpdateCellVectors(Cell.Leaving, AgeLeaving,1);  
    Cell.CellCount = Cell.CellCount + Cell.ToAdd;
    Cell.ToAdd=0;
    Ef.ToAddFromJ=0;
end

if M.ToAdd > 0
   M.CellCount = M.CellCount + M.ToAdd;
   AgeTransit = (steptime.*geornd(pdivM,1,M.ToAdd)) + t + steptime;
   AgeLeaving =  (steptime.*geornd(pleaveM,1,M.ToAdd)) + t + steptime;
   idx = AgeTransit >= AgeLeaving;
   AgeTransit(idx | AgeTransit>stoptime) = [];
   AgeLeaving(~idx | AgeLeaving>stoptime) = [];
   M.Transit = UpdateCellVectors(M.Transit, AgeTransit,1);  
   M.Leaving = UpdateCellVectors(M.Leaving, AgeLeaving,1);  
   M.ToAdd=0;
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