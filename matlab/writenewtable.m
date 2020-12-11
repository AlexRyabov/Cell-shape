function writenewtable(T, filename)
%this function first deletes the old table
delete(filename);
writetable(T, filename);