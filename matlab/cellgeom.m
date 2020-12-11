%% Calculate surface area and volume for cells of various geometry
%% Alexey Ryabov 2020
% the script loads data from '..\data\CellSamples.xlsx'
% and saves results into '..\data\CellSamples_VA.xlsx'
%
% if SpecAndSampleAverage = 1 (see below) then the script first 
% averages linear dimensions for each species and location
% For each record (agregated record in case of averaging)
% the script finds area and volume of cells
% equvivalent spherical radius based on volume (Eqv_Rad_v = (3/4 * V/pi))^(1/3) 
% and equvivalent spherical radius based on surface area (Eqv_Rad_s =
% sqrt(A/(4 * pi))) ;
% surface extension 
% aspect ratio 
% minimal, middle and maximal linear cell dimensions 
% and classify the cell as prolate/oblate or compact

% shape IDs and formulas for surface area and volume 
% are in CalculationsOfCellVolume&Surface.pdf and also 
%  http://phytobioimaging.unisalento.it/Products/AtlasOfShapes.aspx?ID_Tipo=0

clearvars;
tLinDim = readtable('..\data\CellSamples.xlsx');

%get unique shapes 
tUniqShapes = unique(tLinDim(:, {'ShapeCode', 'ShapeType'}), 'rows');


%% this part is a preparation for averaging for each species and each location
SpecAndSampleAverage = 1;  % if 1 then avarage for species and sample, otherwise only for species
if SpecAndSampleAverage
    %averrage for species and sample
    tMeanSizeSpecies = grpstats(tLinDim(:, {...
        'CodeEcosystem', 'Genus_species', 'Genus_species_view_', 'ShapeCode', 'ShapeType', ...
        'HD'    'a'    'b'    'c'    'd'    'h'    'a1'    'b1'    'c1'    'h1'    'd1' ...
        'a2'    'b2'    'c2'    'h2'    'd2'    'a3'    'b3'    'c3'    'h3'    'd3'...
        'h4'    'd4'    'e'    'l'    'l1'    'l2'    'l3'    'l4'    'd11'    'd22'    'd33'    'd44'   }), ...
        {'CodeEcosystem', 'Genus_species', 'Genus_species_view_', 'ShapeCode', 'ShapeType'}, {'mean'});
    for i = 1:size(tMeanSizeSpecies, 2)
        s = tMeanSizeSpecies.Properties.VariableNames{i};
        k = findstr('mean_', s);
        if ~isempty(k) && k == 1
            tMeanSizeSpecies.Properties.VariableNames{i} = s(6:end);
        end
    end
    tLinDim = tMeanSizeSpecies;
end

%Add columns 
tLinDim.Volume = NaN(size(tLinDim, 1), 1);
tLinDim.Area   = NaN(size(tLinDim, 1), 1);
tLinDim.Dim1   = NaN(size(tLinDim, 1), 1);
tLinDim.Dim2   = NaN(size(tLinDim, 1), 1);
tLinDim.Dim3   = NaN(size(tLinDim, 1), 1);


%%
%find surface area, volume, min, med and max dimesions
%using shape code 
for iSh = 1:size(tUniqShapes, 1)
    iSh
    %
    ind = tLinDim.ShapeCode == tUniqShapes.ShapeCode(iSh);
    HD = tLinDim.HD(ind);
    a  = tLinDim.a(ind);
    b  = tLinDim.b(ind);
    c  = tLinDim.c(ind);
    d  = tLinDim.d(ind);
    h  = tLinDim.h(ind);
    a1 = tLinDim.a1(ind);
    b1 = tLinDim.b1(ind);
    c1 = tLinDim.c1(ind);
    h1 = tLinDim.h1(ind);
    d1 = tLinDim.d1(ind);
    a2 = tLinDim.a2(ind);
    b2 = tLinDim.b2(ind);
    c2 = tLinDim.c2(ind);
    h2 = tLinDim.h2(ind);
    d2 = tLinDim.d2(ind);
    a3 = tLinDim.a3(ind);
    b3 = tLinDim.b3(ind);
    c3 = tLinDim.c3(ind);
    h3 = tLinDim.h3(ind);
    d3 = tLinDim.d3(ind);
    h4 = tLinDim.h4(ind);
    d4 = tLinDim.d4(ind);
    e  = tLinDim.e(ind);
    l  = tLinDim.l(ind);
    %%
    switch tUniqShapes.ShapeCode(iSh)
        
        
        case 1 %Sphere
            tLinDim.Volume(ind)= pi/6 * d.^3;
            tLinDim.Area(ind)= pi * d.^2;
            tLinDim.Dim1(ind)= d;
            tLinDim.Dim2(ind)= d;
            tLinDim.Dim3(ind)= d;
        case 2 %Prolate spheroid
            tLinDim.Area(ind)  = pi * d/2 .* (d + h.^2 ./sqrt(h.^2 - d.^2) .* asin(sqrt(h.^2 - d.^2)./h));
            %equation from maple gives exctly the same resutl. TESTED
            %A=(1/2)*pi*d.*(d+2*h.^2.*asin(sqrt(-4*d.^2+4*h.^2)./(2*h))./sqrt(-4*d.^2+4*h.^2))
            tLinDim.Volume(ind)= pi/6 * d.^2 .* h;
            tLinDim.Dim1(ind)=d;
            tLinDim.Dim2(ind)=d;
            tLinDim.Dim3(ind)=h;
        case 3  %Cylinder
            tLinDim.Area(ind)   = pi * d .* (d/2 + h);
            tLinDim.Volume(ind) = pi/4 .* d.^2 .* h;
            tLinDim.Dim1(ind)   = d;
            tLinDim.Dim2(ind)   = d;
            tLinDim.Dim3(ind)   = h;
        case 4 %Ellipsoid
            tLinDim.Area(ind)= pi/4 .*(b+c).*...
                (...
                (b+c)/2 + 2 * h.^2./sqrt(4 * h.^2 - (b + c).^2) .* asin(sqrt(4 * h.^2 - (b + c).^2)./(2 * h) ) ...
                );
            tLinDim.Volume(ind)  = pi/6 * b.*c.*h;
            tLinDim.Dim1(ind)=b;
            tLinDim.Dim2(ind)=c;
            tLinDim.Dim3(ind)=h;
        case 5 %Cone
            tLinDim.Area(ind)   = pi/4 * d .* (d + sqrt(4 * h.^2 + d.^2));
            tLinDim.Volume(ind) = pi/12 * d.^2 .*h;
            tLinDim.Dim1(ind)   = d;
            tLinDim.Dim2(ind)   = d;
            tLinDim.Dim3(ind)   = h;
        case 7 %Parallelepiped
            tLinDim.Area(ind)   = 2*a.*b + 2 * b.*c + 2 * a.*c;
            tLinDim.Volume(ind) = a.*b.*c;
            tLinDim.Dim1(ind)   = a;
            tLinDim.Dim2(ind)   = b;
            tLinDim.Dim3(ind)   = c;
        case 8 %Prism on elliptic base
            %tLinDim.Area(ind)   = pi/2 * (a.*b + c .* (a + b));
            %this formula from http://phytobioimaging.unisalento.it/en-us/products/AtlasOfShapes.aspx?ID_Tipo=3
            %suggests a first order approximation for the Area, better this
            tLinDim.Area(ind)   =c*pi.*((1/2)*a+(1/2)*b).*(1+(a-b).^2./(4*(a+b).^2))+(1/2)*pi*a.*b;
            tLinDim.Volume(ind) = (1/4)*pi*a.*b.*c;
            tLinDim.Dim1(ind)   = a;
            tLinDim.Dim2(ind)   = b;
            tLinDim.Dim3(ind)   = c;
        case 9 %Prism on parallelogram base
            %tLinDim.Area(ind)   = a.*b + sqrt(a.^2 + b.^2)/4 .* c; %this
            %formula from Hillbrand 1999 is wrong, it gives approximately 50% error
            tLinDim.Area(ind)   = a.*b+2*sqrt(a.^2+b.^2).*c;
            tLinDim.Volume(ind) = (1/2)*a.*b.*c;
            tLinDim.Dim1(ind)   = a;
            tLinDim.Dim2(ind)   = b;
            tLinDim.Dim3(ind)   = c;
        case 10  %Cube
            tLinDim.Area(ind)   = 6 * a.^2;
            tLinDim.Volume(ind) = a.^3;
            tLinDim.Dim1(ind)   = a;
            tLinDim.Dim2(ind)   = a;
            tLinDim.Dim3(ind)   = a;
        case 11  %Prism on triangle base 1
            %there is a typo in the formula for area 3 ab + \sqrt(3)/2 *
            %a^2 on web, it should be
            tLinDim.Area(ind)   = 3.*a.*c+(1/2)*sqrt(3).*a.^2;
            tLinDim.Volume(ind) = (1/4)*sqrt(3)*a.^2.*c;
            tLinDim.Dim1(ind)   = a;
            tLinDim.Dim2(ind)   = sqrt(3)/2*a;  %height
            tLinDim.Dim3(ind)   = c;
        case 12  %Half prism on elliptic base
            % A = pi/4 * (a .* b + b .* c  + a.*c ) + a.*c;
            %this formula from http://phytobioimaging.unisalento.it/en-us/products/AtlasOfShapes.aspx?ID_Tipo=3
            %is wrong (also in Hillebrand 1999)
            % one can use A = pi/2 * a.*b + pi/2 * (a/2 + b).*c + a.*c;
            % a more presise formula is
            tLinDim.Area(ind)   = (1/2)*pi*((1/2)*a+b).*(1+(a-2*b).^2./(4*(a+2*b).^2)).*c+a.*c+(1/2)*pi*a.*b;
            tLinDim.Volume(ind) = (1/4)*pi*a.*b.*c;
            tLinDim.Dim1(ind)   = a;
            tLinDim.Dim2(ind)   = b;
            tLinDim.Dim3(ind)   = c;
        case 14  %Double cone
            tLinDim.Area(ind)   = (1/2)*pi*d.*sqrt(d.^2+h.^2);
            tLinDim.Volume(ind) = (1/12)*pi*d.^2.*h;
            tLinDim.Dim1(ind)   = d;
            tLinDim.Dim2(ind)   = d;
            tLinDim.Dim3(ind)   = h;
        case 15 %2 truncated cones. I didn't check it (there is only one cell)
            l = sqrt((d2/2-d1/2).^2 + h.^2);
            tLinDim.Area(ind)   = pi * l .* (d1 + d2) + pi/2 * d1.^2;
            tLinDim.Volume(ind) = pi/6 * h .* (d1.^2 + d1.*d2 + d2.^2);
            tLinDim.Dim1(ind)   = max(d1, d2);
            tLinDim.Dim2(ind)   = max(d1, d2);
            tLinDim.Dim3(ind)   = 2*h;
        case 16 %Prolate spheroid + 2 Cylinders
            tLinDim.Area(ind)   = pi*d1.*h1+pi*d2.*h2+pi*d.*h.^2.*asin(sqrt(-d.^2+h.^2)./h)./(2*sqrt(-d.^2+h.^2))+(1/2)*pi*d.^2;
            tLinDim.Volume(ind) = (1/4)*pi*d1.^2.*h1+(1/4)*pi*d2.^2.*h2+(1/6)*pi*d.^2.*h;
            tLinDim.Dim1(ind)   = d;
            tLinDim.Dim2(ind)   = d;
            tLinDim.Dim3(ind)   = h1 + h2 + h;
        case 17  %Cylinder +2 cones % I assume that the cones are equilateral
            %Both Hillbrand 1999 and web have mistakes in the formulas
            %tLinDim.Area(ind)   = pi*d.^2+pi*d.*(h-sqrt(3)*d);
            %tLinDim.Volume(ind) = (1/12)*pi*d.^3*sqrt(3)+(1/4)*pi*d.^2.*(h-sqrt(3)*d);
            %here I assume that the cone height equals d/2
            tLinDim.Area(ind)   = (1/2)*pi*d.^2*sqrt(2)+pi*d.*(h-d);
            tLinDim.Volume(ind) = -(1/6)*pi*d.^3+(1/4)*pi*d.^2.*h;
            tLinDim.Dim1(ind)   = d;
            tLinDim.Dim2(ind)   = d;
            tLinDim.Dim3(ind)   = h;
        case 19 % cone+half sphere. The Volume is wrong on web, Area is wrong on web
            tLinDim.Area(ind)   = (1/4)*pi*d.*(sqrt(2*d.^2-4*d.*h+4*h.^2)+2*d);
            tLinDim.Volume(ind) = (1/24)*pi*d.^2.*(d+2*h);
            tLinDim.Dim1(ind)   = d;
            tLinDim.Dim2(ind)   = d;
            tLinDim.Dim3(ind)   = h;
        case 20  %Half ellipsoid + Cone (on elliptic base)
            %on web there was only Volume
            ConeSideArea = (1/2)*pi*((1/4)*b.*sqrt(b.^2+4*h1.^2)+(1/4)*c.*sqrt(c.^2+4*h1.^2));
            HalfEllArea = (1/8)*pi*(b+c).*(...
                (1/2)*b+(1/2)*c+2*h.^2.*asin(sqrt(4*h.^2-(b+c).^2)./(2*h))./sqrt(4*h.^2-(b+c).^2)...
                );
            tLinDim.Area(ind)   = ConeSideArea + HalfEllArea;
            tLinDim.Volume(ind) = (1/12)*pi*c.*b.*h1+(1/12)*pi.*b.*c.*h;
            tLinDim.Dim1(ind)   = b;
            tLinDim.Dim2(ind)   = c;
            tLinDim.Dim3(ind)   = h/2 + h1;
        case 21 %Prism on elliptic base+ box, ASh gives a correct formula for V, but a bit wrong for A
            P = pi*((1/2)*a+(1/2)*b).*(1+(a-b).^2./(4*(a+b).^2))+2*a1;
            tLinDim.Area(ind)   =P.*c +2*a1.*b1+(1/2)*pi*a.*b;
            tLinDim.Volume(ind) = a1.*b1.*c+(1/4)*pi*a.*b.*c;
            tLinDim.Dim1(ind)   = b;
            tLinDim.Dim2(ind)   = c;
            tLinDim.Dim3(ind)   = a+a1;
        case 22 %Cylinder + 2 Half spheres
            tLinDim.Area(ind)   = d*pi.*(d+h);
            tLinDim.Volume(ind) = (1/6)*pi*d.^3+(1/4)*pi*d.^2.*h;
            tLinDim.Dim1(ind)   = d;
            tLinDim.Dim2(ind)   = d;
            tLinDim.Dim3(ind)   = h+d;
        case 23 %Ellipsoid+2cones+cylinder
            tLinDim.Area(ind)   = (1/4)*pi*(b+c).*(...
                (1/2)*b+(1/2)*c+2*h.^2.*asin(sqrt(4*h.^2-(b+c).^2)./(2*h))./sqrt(4*h.^2-(b+c).^2)...
                ) ...
                -(1/4)*pi*d2.^2-(1/4)*pi*d3.^2+h1.*d1*pi+(1/2)*pi.*d2.*sqrt(h2.^2+(1/4)*d2.^2)...
                +(1/2)*pi.*d3.*sqrt(h3.^2+(1/4)*d3.^2);
            tLinDim.Volume(ind) = (1/6)*pi*b.*c.*h+(1/4)*pi*d1.^2.*h1+(1/12)*pi*d2.^2.*h2+(1/12)*pi*d3.^2.*h3;
            tLinDim.Dim1(ind)   = b;
            tLinDim.Dim2(ind)   = c;
            tLinDim.Dim3(ind)   = h+h1+max(h2,h3);
        case 24 %Ellipsoid + Cone
            tLinDim.Area(ind)   = ...
                (1/4)*pi*(b+c).*(...
                (1/2)*b+(1/2)*c+2*h.^2.*asin(sqrt(4*h.^2-(b+c).^2)./(2*h))./sqrt(4*h.^2-(b+c).^2)...
                )...
                -(1/4)*pi*d1.^2+(1/2)*pi*d1.*sqrt(h1.^2+(1/4)*d1.^2);
            tLinDim.Volume(ind) = (1/6)*pi*b.*c.*h+(1/12)*pi*d1.^2.*h1;
            tLinDim.Dim1(ind)   = b;
            tLinDim.Dim2(ind)   = c;
            tLinDim.Dim3(ind)   = h+h1;
        case 25 %Cylinder + 3 Cones
            tLinDim.Area(ind)   = pi*((1/2)*d1+(1/2)*d4).*sqrt(((1/2)*d4-(1/2)*d1).^2+h4.^2)...
                +(1/4)*pi*d4.^2+h1.*d1*pi+(1/4)*pi*d1.^2-(1/4)*pi*d3.^2-(1/4)*pi*d2.^2 ...
                +(1/2)*pi*d3.*sqrt(h3.^2+(1/4)*d3.^2)+(1/2)*pi*d2.*sqrt(h2.^2+(1/4)*d2.^2);
            tLinDim.Volume(ind) = (1/12)*pi*h4.*(d1.^2+d1.*d4+d4.^2)...
                +(1/4)*pi*d1.^2.*h1+(1/12)*pi*d2.^2.*h2+(1/12)*pi*d3.^2.*h3;
            tLinDim.Dim1(ind)   = d4;
            tLinDim.Dim2(ind)   = d4;
            tLinDim.Dim3(ind)   = h1+h4+max(h2, h3);
        case 27 %Half sphere
            tLinDim.Area(ind)   = 3*pi*d.^2*(1/4);
            tLinDim.Volume(ind) = (1/12)*pi*d.^3;
            tLinDim.Dim1(ind)   = d;
            tLinDim.Dim2(ind)   = d;
            tLinDim.Dim3(ind)   = d/2;
        case 34 %2 half ellipsoids + priism on elliptic base
            tLinDim.Area(ind)   = ...
                (1/8)*pi*(2*b1+c1).*(...
                b1+(1/2)*c1+2*h1.^2.*asin(sqrt(4*h1.^2-(2*b1+c1).^2)./(2*h1))./sqrt(4*h1.^2-(2*b1+c1).^2)...
                )...
                +(1/4)*pi*h1.*c1-(1/4)*pi*a.*c1...
                +(1/8)*pi*(2*b2+c2).*(...
                b2+(1/2)*c2+2*h2.^2.*asin(sqrt(4*h2.^2-(2*b2+c2).^2)./(2*h2))./sqrt(4*h2.^2-(2*b2+c2).^2)...
                )...
                +(1/4)*pi*h2.*c2-(1/4)*pi*a.*c2+pi*((1/2)*a+(1/2)*c1).*(1+(a-c1).^2./(4*(a+c1).^2)).*c;            tLinDim.Volume(ind) = (1/12)*pi*b1.*c1.*h1+(1/12)*pi*b2.*c2.*h2+(1/4)*pi*a.*b.*c;
            tLinDim.Volume(ind)   = (1/6)*pi*b1.*c1.*h1+(1/6)*pi*b2.*c2.*h2+(1/4)*pi*a.*c1.*c;
            tLinDim.Dim1(ind)   = mean([c1, c2], 2);
            tLinDim.Dim2(ind)   = b1+c+b2;
            tLinDim.Dim3(ind)   = mean([h1, h2], 2);
        case 35  %Cymbelloid. in Hi99 we dont have area, in ASh the area is wrongly found
            %(it should be b instead of c and arcsin(beta))
            tLinDim.Area(ind)   = b.*(2*b+2*a.^2.*asin(sqrt(4*a.^2-16*b.^2)./(2*a))./sqrt(4*a.^2-16*b.^2))...
                .*asin(c./(2*b))+(1/2)*pi*a.*b;
            tLinDim.Volume(ind) = 2*a.*b.^2.*asin(c./(2*b))*(1/3);
            tLinDim.Dim1(ind)   = b;
            tLinDim.Dim2(ind)   = c;
            tLinDim.Dim3(ind)   = a;
        case 40 %Gomphonemoid
            tLinDim.Area(ind)   = (1/2)*b.*(2*h+pi*h.*asin(c./(2*h))+((1/2)*pi-2)*b);
            tLinDim.Volume(ind) = (1/4)*h.*b.*(h+((1/4)*pi-1)*b).*asin(c./(2*h));
            tLinDim.Dim1(ind)   = b;
            tLinDim.Dim2(ind)   = c;
            tLinDim.Dim3(ind)   = h;
        case 41  %Sickle-shaped prism
            b = a;
            b2 = 0 * b; %%assume the inner semi axis equals 0
            tLinDim.Area(ind)   = (1/2)*pi*(b.*c+b.*h+b2.*c-b2.*h+c.*h);
            tLinDim.Volume(ind) = (1/4)*pi*c.*h.*(a);
            tLinDim.Dim1(ind)   = a;
            tLinDim.Dim2(ind)   = c;
            tLinDim.Dim3(ind)   = h;
        case 43  %Prism on elliptic base + 4 Cones
            tLinDim.Area(ind)  = c.*pi.*((1/2)*a+(1/2)*b).*(1+(a-b).^2./(4*(a+b).^2)) ...
                +(1/2)*pi*a.*b-(1/4)*pi*d1.^2-(1/4)*pi*d2.^2-(1/4)*pi*d3.^2-(1/4)*pi*d4.^2 ...
                +(1/2)*pi*d1.*sqrt(h1.^2+(1/4)*d1.^2)+(1/2)*pi*d2.*sqrt(h2.^2+(1/4)*d2.^2) ...
                +(1/2)*pi*d3.*sqrt(h3.^2+(1/4)*d3.^2)+(1/2)*pi*d4.*sqrt(h4.^2+(1/4)*d4.^2);
            tLinDim.Volume(ind) = (1/4)*pi*a.*b.*c+(1/12)*pi*d1.^2.*h1+(1/12)*pi*d2.^2.*h2 ...
                +(1/12)*pi*d3.^2.*h3+(1/12)*pi*d4.^2.*h4;
            tLinDim.Dim1(ind)   = a;
            tLinDim.Dim2(ind)   = b;
            tLinDim.Dim3(ind)   = max(h1, h2) + c + max(h3, h4);
        case 44  %Pyramid (rectangular base)
            tLinDim.Area(ind)   = sqrt(d.^2+4*h.^2).*d+d.^2;
            tLinDim.Volume(ind) = (1/3)*d.^2.*h;
            tLinDim.Dim1(ind)   = d;
            tLinDim.Dim2(ind)   = d;
            tLinDim.Dim3(ind)   = h;
        case 46  %Prisma on triangle-base 2
            tLinDim.Area(ind)   = 3*a.*b+(1/2)*a.^2*sqrt(3);
            tLinDim.Volume(ind) = (1/4)*a.^2.*sqrt(3).*b;
            tLinDim.Dim1(ind)   = a;
            tLinDim.Dim2(ind)   = sqrt(3)/2 * a;
            tLinDim.Dim3(ind)   = b;
        case 51 %2 Half ellipsoids
            tLinDim.Area(ind)   = (1/8)*pi*(2*b1+c1).*(...
                b1+(1/2)*c1+2*h1.^2.*asin(sqrt(4*h1.^2-(2*b1+c1).^2)./(2*h1))./sqrt(4*h1.^2-(2*b1+c1).^2)...
                )...
                +(1/8)*pi*(2*b2+c2).*(...
                b2+(1/2)*c2+2*h2.^2.*asin(sqrt(4*h2.^2-(2*b2+c2).^2)./(2*h2))./sqrt(4*h2.^2-(2*b2+c2).^2));
            tLinDim.Volume(ind) = (1/6)*pi*b1.*c1.*h1+(1/6)*pi*b2.*c2.*h2;
            tLinDim.Dim1(ind)   = mean([c1, c2], 2);
            tLinDim.Dim2(ind)   = b1+b2;
            tLinDim.Dim3(ind)   = mean([h1, h2], 2);
        otherwise
            error('Unknown shape ID')
    end
    
end

%Surface extensions (inverse sphericity)
tLinDim.Eqv_Rad_v = (tLinDim.Volume/(4/3 * pi)).^(1/3);
tLinDim.Eqv_Rad_s = (tLinDim.Area/(4 * pi)).^(1/2);
tLinDim.Surf_ext  = tLinDim.Area./(4 * pi * tLinDim.Eqv_Rad_v.^2);

%Add mininal, middle and maximal dimensions
DimsSorted = [tLinDim.Dim1, tLinDim.Dim2, tLinDim.Dim3];
DimsSorted = [sort(DimsSorted, 2), DimsSorted];
tLinDim.LMin = DimsSorted(:, 1);
tLinDim.LMid = DimsSorted(:, 2);
tLinDim.LMax = DimsSorted(:, 3);

%ratio of the middle dimension to the minimal and maximal
tLinDim.DmMid2Min = tLinDim.LMid./tLinDim.LMin;
tLinDim.DmMax2Mid = tLinDim.LMax./tLinDim.LMid;


%Acspect ratio
tLinDim.AspectRatio = tLinDim.LMax./tLinDim.LMin;  %Maximal 
ind = tLinDim.DmMax2Mid<tLinDim.DmMid2Min;
tLinDim.AspectRatio(ind) = 1./tLinDim.AspectRatio(ind);
%Prolate/oblate
tLinDim.Prolate = ~ind;

%Prolate/Oblate/compact
CompactCellRatio = 1.5;
ind = tLinDim.DmMax2Mid<CompactCellRatio & tLinDim.DmMid2Min<CompactCellRatio;
tLinDim.Elongation = repmat(['U'], size(tLinDim, 1), 1);  %unkknown cells
tLinDim.Elongation(ind) = 'C';                              %compact cells
tLinDim.Elongation(~ind & (tLinDim.DmMax2Mid>=tLinDim.DmMid2Min)) = 'P';  %prolate
tLinDim.Elongation(~ind & (tLinDim.DmMax2Mid<tLinDim.DmMid2Min)) = 'O';   %oblate







%%
%save data
writenewtable(tLinDim,'..\data\CellSamples_VA.xlsx')
