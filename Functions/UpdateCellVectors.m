function [updated_warehouse] = UpdateCellVectors(inventory, shipment,sign)

updated_warehouse=inventory;

if ~isempty(shipment)
    if isempty(shipment)
        updated_warehouse = inventory;
        return;
    end
    shipment = sort(shipment);
    diffs = [1, diff(shipment) ~= 0];
    unique_ids = shipment(diffs == 1);
    counts = diff([find(diffs), length(shipment) + 1]) * sign;
    unique_ids = reshape(unique_ids, 1, []);
    counts = reshape(counts, 1, []);

    if isempty(inventory)
        if sign == 1
            updated_warehouse = [unique_ids; max(counts, 0)];
        else
            updated_warehouse = [];
        end
        return;
    end
    all_ids = [inventory(1, :) unique_ids];
    all_counts = [inventory(2, :) counts];
    [sorted_ids, sort_idx] = sort(all_ids);
    sorted_counts = all_counts(sort_idx);
    diffs = [1, diff(sorted_ids) ~= 0];
    final_ids = sorted_ids(diffs == 1);
    final_counts = accumarray(cumsum(diffs)', sorted_counts);
    valid_idx = final_counts > 0;
    updated_warehouse = [final_ids(valid_idx); final_counts(valid_idx)'];
end

end

