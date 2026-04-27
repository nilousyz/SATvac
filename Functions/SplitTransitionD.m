function [CellIn,CellOut] = SplitTransitionD(CellIn,CellOut,pSplit,CellLeavingProb,t,steptime,stoptime)

if ~isempty(CellIn.Transit)
    idx = CellIn.Transit(1,:) < t + 0.1;
    CellIn.CellCount = CellIn.CellCount - sum(CellIn.Transit(2,idx));   
    CellOut.ToAdd = CellOut.ToAdd + (2*sum(CellIn.Transit(2,idx)));     
    CellIn.Transit = CellIn.Transit(:,~idx);                            

end

if  CellOut.ToAdd > 0
    AgeLeaving = (steptime.*geoinv(rand(1,CellOut.ToAdd),CellLeavingProb) + t + steptime);
    AgeTransit = (inverseCDF(CellOut.ToAdd,pSplit) + t + steptime);
    idx = AgeTransit >= AgeLeaving;
    AgeTransit(idx | AgeTransit>stoptime) = [];
    AgeLeaving(~idx | AgeLeaving>stoptime) = [];
    CellOut.Transit = UpdateCellVectors(CellOut.Transit, AgeTransit,1);  
    CellOut.Leaving = UpdateCellVectors(CellOut.Leaving, AgeLeaving,1);  
    CellOut.CellCount = CellOut.CellCount + CellOut.ToAdd;
    CellOut.ToAdd=0;
end
end





