%%Alex Ryabov, 2019, 
%%replaces unbreacable spaces with normal
%%Remove trailing spaces at the end and 
%% tTable.StringColumn = f_TrimColumn(tTable.StringColumn);

function col = f_TrimColumn(col)
for i = 1:numel(col)

    i
    a = col{i};
    ind = double(a) == 160;
    if sum(ind)
       a(ind) = ' ';
    end
    sp_ind = isspace(a);
    while sp_ind(end)
        a = a(1:end-1);
        sp_ind = isspace(a);
    end
    col{i} = a;
end
