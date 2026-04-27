function [Cell,K] = CellLeavingToNewK(Cell,K,t,steptime)

if  ~isempty(Cell.Leaving)
    idx = Cell.Leaving(1,:) < t + 0.1;
    Cell.CellCount = Cell.CellCount - sum(Cell.Leaving(2,idx));
    K.ToAdd = sum(Cell.Leaving(2,idx));
    K.CellCount = K.CellCount + K.ToAdd;
    Cell.Leaving = Cell.Leaving(:,~idx);

    if K.ToAdd > 50000000
        warning('Excessive K')
        error('Error Occurred')
    end
    if  K.ToAdd > 0
        K.ToAdd = 0;
    end

end

end

