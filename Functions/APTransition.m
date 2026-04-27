function [CellIn,CellOut] = APTransition(CellIn,CellOut,pBP,t,steptime)

if ~isempty(CellIn.Transit)
    idx = CellIn.Transit(1,:) < t+0.1;
    CellIn.CellCount = CellIn.CellCount - sum(CellIn.Transit(2,idx));
    CellOut.ToAdd = CellOut.ToAdd + sum(CellIn.Transit(2,idx));
    CellIn.Transit = CellIn.Transit(:,~idx);
end

if  CellOut.ToAdd > 0
    AgeTransit = (inverseCDF(CellOut.ToAdd,pBP) + t + steptime);
    CellOut.Transit = UpdateCellVectors(CellOut.Transit, AgeTransit,1);
    CellOut.CellCount = CellOut.CellCount + CellOut.ToAdd;
    CellOut.ToAdd=0;
end
end


