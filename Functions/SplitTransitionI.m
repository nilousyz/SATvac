function [CellIn,CellOut] = SplitTransitionI(CellIn,CellOut,t)

if ~isempty(CellIn.Transit)
    idx = CellIn.Transit(1,:) < t+0.1;
    CellIn.CellCount = CellIn.CellCount - sum(CellIn.Transit(2,idx));       
    CellOut.ToAdd = CellOut.ToAdd + (2*sum(CellIn.Transit(2,idx)));         
    CellIn.Transit = CellIn.Transit(:,~idx);                                
    
end



end