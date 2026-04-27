function [CellIn,CellOut] = SplitTransitionC(CellIn,CellOut,pSplit,t,steptime,stoptime)
if ~isempty(CellIn.Transit)
    idx = CellIn.Transit(1,:) < t + 0.1;
    CellIn.CellCount = CellIn.CellCount - sum(CellIn.Transit(2,idx));       
    CellOut.ToAdd = CellOut.ToAdd + (2*sum(CellIn.Transit(2,idx)));         
    CellIn.Transit = CellIn.Transit(:,~idx);                                
end

if  CellOut.ToAdd > 0
    AgeTransit = (inverseCDF(CellOut.ToAdd,pSplit) + t + steptime);
    idx = AgeTransit > stoptime;
    AgeTransit(idx) = [];
    CellOut.Transit = UpdateCellVectors(CellOut.Transit, AgeTransit,1);  
    CellOut.CellCount = CellOut.CellCount + CellOut.ToAdd;
    CellOut.ToAdd=0;
end
end





