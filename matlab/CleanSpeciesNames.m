function tTable = CleanSpeciesNames(tTable, AllowedLostRows)
%% Alexey Ryabov, 2020
%%changes genus names using Genus.Species_lea.xlsx, where the names of
%%genus were taken from AlgaeDatabase
%if AllowedLostRows > 0 then some rows might be not found in the data on
%abundance

if ~exist('AllowedLostRows', 'var')
    AllowedLostRows = 0;
end

%%correct species names
%remove last space from species names
tTable.SpeciesName = f_TrimColumn(tTable.SpeciesName);

%%adds species names from Genus.Species_lea.xlsx
tGenusSpeciesMatch = readtable('..\data\Genus.Species_lea.xlsx');
tGenusSpeciesMatch.OurGenusName = f_TrimColumn(tGenusSpeciesMatch.OurGenusName);
tTable1 = innerjoin(tTable, tGenusSpeciesMatch(:, {'OurGenusName', 'AlgaeDataBaseName'}), ...
    'LeftKeys', 'SpeciesName', 'RightKey', 'OurGenusName');
if (size(tTable1, 1) > size(tTable, 1)) ||  (size(tTable1, 1) < size(tTable, 1) - AllowedLostRows)
    error('Some data were lost');
end
tTable1.SpeciesName = tTable1.AlgaeDataBaseName;
tTable1.AlgaeDataBaseName = [];
tTable = tTable1;

end

