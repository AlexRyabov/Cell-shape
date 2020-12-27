"""
receives the fieldnames and dimension values for a single species
calculates parameters like Volume, Area, and returns them together with sorted Dimensions (min, mid, max)
Disclaimer: calculations of A and V are pythonized from the matlab script cellgeom.m by A. Ryabov
Tested for: python 3.6
"""

# Library imports
import numpy as np
import re

__author__ = "Onur Kerimoglu"
__email__ = "onur.kerimoglu@uol.de"
__credits__ = ["Onur Kerimoglu", "Alexey Ryabov"]
__version__ = "1.0.0"  # December 2020

# global parameters and aliases
pi = np.pi
sqrt = np.sqrt
asin = np.arcsin
nan = np.nan

def calc_geom(fields,headersin,dataset):
    shapecol = headersin.index('Geometric_shape')
    gshape = fields[shapecol]
    if gshape[0] == ' ':
        gshape = gshape[1:]
    gshape.replace('  ', '_')
    gshape.replace('   ', '_')
    #todo: decapitalize
    Found = True #default: assume the shape is recognized and dimensions can be extracted
    if (gshape in ['1','sphere','sphere-10%','sphere-20%','sphere-25%']):
        # In Olenina dataset: sphere/-10%/-20%/-25% (regex search does not work because of confounding instances)
        if dataset=='Olenina':
            d = get_dim(fields, headersin, 'D1',dataset)
        elif dataset=='Roselli':
            Found = False
            unable2extrdim(dataset,gshape)
        if Found:
            V,A,Dsorted = shape_1(d)
    elif (gshape in ['2','prolate spheroid']) or (re.search('rotational ellipsoid',gshape)):
        # In Olenina dataset (=prolate spheroid): rotational ellipsoid/ x 0.5/-20%
        if dataset == 'Olenina':
            d = get_dim(fields, headersin, 'D1',dataset) #check: /2 ?
            h = get_dim(fields, headersin, 'L1',dataset) #check: /2 ?
            if h<d: #this results in sqrt of negative number for area calculation (sqrt (h^2 - d^2))
                #solution: swap dimensions (check)
                d=get_dim(fields, headersin, 'L1',dataset)
                h=get_dim(fields, headersin, 'D1',dataset)
        elif dataset=='Roselli':
            Found = False
            unable2extrdim(dataset,gshape)
        if Found:
            if d == h:
                gshape='sphere'
                V, A, Dsorted = shape_1(d) #calculating for shape_2 with d==h results in division by 0
            else:
                V, A, Dsorted = shape_2(d,h)
    elif gshape in ['3', 'cylinder']:
        if dataset == 'Olenina':
            d = get_dim(fields, headersin, 'D1',dataset)
            h = get_dim(fields, headersin, 'H',dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset,gshape)
        if Found:
            V, A, Dsorted = shape_3(d, h)
    elif (gshape in ['4', 'ellipsoid']) or (re.search('flattened ellipsoid',gshape) is not None):
        # In Olenina dataset (=ellipsoid): flattened ellipsoid/ - 20%/-20%
        if dataset == 'Olenina':
            b = get_dim(fields, headersin, 'D1') #check: / 2
            c = get_dim(fields, headersin, 'D2') #check: / 2
            h = get_dim(fields, headersin, 'L1') #check: / 2
            if np.isnan(h):
                h = get_dim(fields, headersin, 'H')  # check: / 2
            if np.isnan(h):
                print('!')
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset,gshape)
        if Found:
            V, A, Dsorted = shape_4(b, c, h)
    elif (gshape in ['5', 'cone', 'cone-10%']):
        # In Olenina dataset: cone/-10%
        if dataset == 'Olenina':
            d = get_dim(fields, headersin, 'D1')
            h = get_dim(fields, headersin, 'H')
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset,gshape)
        if Found:
            V, A, Dsorted = shape_5(d, h)
    elif (gshape in ['7','parallelepiped']) or (re.search('parallelipiped-',gshape) is not None):
        # In Olenina dataset: parallelipiped/-10%/-20%/-25%/-30%/-40%
        if dataset == 'Olenina':
            a = get_dim(fields, headersin, 'H')
            b = get_dim(fields, headersin, 'L1')
            c = get_dim(fields, headersin, 'W')
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset,gshape)
        if Found:
            V, A, Dsorted = shape_7(a, b, c)
    elif (gshape in ['8','prism on elliptic base']) or (re.search('oval cylinder', gshape)):
        # In Olenina dataset (=prism on elliptic base): oval cylinder/-30%
        if dataset == 'Olenina':
            a = get_dim(fields, headersin, 'D1')
            b = get_dim(fields, headersin, 'D2')
            c = get_dim(fields, headersin, 'H')
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset,gshape)
        if Found:
            V, A, Dsorted = shape_8(a, b, c)
    elif gshape in ['9','prism on parallelogram base'] :
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape, dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_9(a, b, c)
    elif gshape in ['10','cube']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape, dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_10(a)
    elif gshape in ['11','prism on triangle base 1', 'parallelipiped/2', 'half parallelipiped']:
        # In Olenina dataset: parallelipiped/2, half parallelipiped
        if dataset == 'Olenina':
            #Check: for these shapes, 3 parameters are provided in Olenina datasets
            H = get_dim(fields, headersin, 'H')
            L1 = get_dim(fields, headersin, 'L1')
            #W = get_dim(fields, headersin, 'W')
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset,gshape)
        if Found:
            V, A, Dsorted = shape_11(a, c)
    elif gshape in ['12','half prism on elliptic base']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape, dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_12(a, b, c)
    elif gshape in ['14','double cone', '2 cones-30%', '2 cones']:
        # In Olenina dataset: 2 cones-30%, 2 cones
        if dataset == 'Olenina':
            d = get_dim(fields, headersin, 'D1')
            h = get_dim(fields, headersin, 'H')
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset,gshape)
        if Found:
            V, A, Dsorted = shape_14(d, h)
    elif gshape in ['15','2 truncated cones', 'two truncated cones']:
        if dataset == 'Olenina':
            Found = False
            unable2extrdim(dataset, gshape)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset,gshape)
        if Found:
            V, A, Dsorted = shape_15(d1, d2, h)
    elif gshape in ['16','prolate spheroid + 2 Cylinders']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape,dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_16(d,d1,d2,h,h1,h2)
    elif gshape in ['17','cylinder + 2 cones']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape,dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_17(d,h)
    elif (gshape in ['19', 'cone + half sphere']) or (re.search('cone + half sphere', gshape)):
        # In Olenina dataset: cone + half sphere/-20%/25%/40%
        if dataset == 'Olenina':
            d = get_dim(fields, headersin, 'D1')
            h = get_dim(fields, headersin, 'H')
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset,gshape)
        if Found:
            V, A, Dsorted = shape_19(d,h)
    elif gshape in ['20', 'half ellipsoid + cone']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape,dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_20(b,c,h,h1)
    elif gshape in ['21','prism on elliptic base+ box']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape,dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_21(a,a1,b,b1,c)
    elif gshape in ['22', 'cylinder + 2 half spheres']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape,dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_22(d,h)
    elif gshape in ['23', 'ellipsoid+2cones+cylinder']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape,dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_23(b, c, d1, d2, d3, h, h1, h2, h3)
    elif gshape in ['24', 'ellipsoid + cone']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape,dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_24(b,c,d1,h,h1)
    elif gshape in ['25', 'cylinder + 3 cones']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape,dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_25(d1,d2,d3,d4,h1,h2,h3,h4)
    elif gshape in ['27', 'half sphere', 'half sphere-30%']:
        if dataset == 'Olenina':
            #Check: For half sphere, 2 parameters are provided in the Olenina dataset: D1 and H
            d = get_dim(fields, headersin, 'D1')
            #H = get_dim(fields, headersin, 'H')
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_27(d)
    elif gshape in ['34', '2 half ellipsoids + prism on elliptic base']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape,dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_34(a,b1,b2,c,c1,c2,h1,h2)
    elif gshape in ['35', 'cymbelloid']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape,dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_35(a,b,c)
    elif gshape in ['38', 'half cone']:
        if dataset == 'Olenina':
            d = get_dim(fields, headersin, 'D1')
            h = get_dim(fields, headersin, 'H')
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_38(d,h) # WARNING: Was not available in cellgeom.m!
    elif gshape in ['40', 'gomphonemoid']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape,dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_40(b,c,h)
    elif gshape in ['41', 'sickle-shaped prism']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape,dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_41(a,c,h)
    elif gshape in ['43', 'prism on elliptic base + 4 cones']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape,dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset,gshape)
        if Found:
            V, A, Dsorted = shape_43(a,b,c,d1,d2,d3,d4,h1,h2,h3,h4)
    elif gshape in ['44', 'pyramid']:
        if dataset == 'Olenina':
            Found = False
            shapenotrec(gshape,dataset)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_44(d,h)
    elif gshape in ['46', 'prism on triangular base', 'prism on triangle base 2',]:
        if dataset == 'Olenina':
            b = get_dim(fields, headersin, 'H')
            a = get_dim(fields, headersin, 'L1')  # in almost every case W=L (Hillebrand99)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_46(a,b)
    elif gshape in ['51', '2 half ellipsoids', '2 rotational ellipsoids']:
        if dataset == 'Olenina':
            Found = False
            unable2extrdim(dataset, gshape)
        elif dataset == 'Roselli':
            Found = False
            unable2extrdim(dataset, gshape)
        if Found:
            V, A, Dsorted = shape_51(b1,b2,c1,c2,h1,h2)
    else:
        Found = False

    if not Found:
        Warning('Unrecognized shape: '+gshape)
        V = nan; A = nan; Dsorted = [nan,nan,nan]

    return gshape, V, A, Dsorted


def unable2extrdim(dataset,gshape):
    Warning('Extracting dims for %s from %s datasets not yet implemented' %(gshape,dataset))


def shapenotrec(dataset,gshape):
    raise(ValueError('%s is not recognized for %s dataset' %(gshape,dataset)))


def shape_1(d):  # Sphere
    V = pi / 6 * d ** 3
    A = pi * d ** 2
    D1 = d
    D2 = d
    D3 = d
    return V, A, np.sort([D1, D2, D3])


def shape_2(d, h):  # Prolate spheroid
    A = pi * d / 2 * (d + h ** 2 / sqrt(h ** 2 - d ** 2) * asin(sqrt(h ** 2 - d ** 2) / h))
    # equation from maple gives exactly the same result. TESTED
    # A=(1/2)*pi*d*(d+2*h**2*asin(sqrt(-4*d**2+4*h**2)/(2*h))/sqrt(-4*d**2+4*h**2))
    V = pi / 6 * d ** 2 * h
    D1 = d
    D2 = d
    D3 = h
    return V, A, np.sort([D1, D2, D3])


def shape_3(d, h):  # Cylinder
    A = pi * d * (d / 2 + h)
    V = pi / 4 * d ** 2 * h
    D1 = d
    D2 = d
    D3 = h
    return V, A, np.sort([D1, D2, D3])


def shape_4(b, c, h):  # Ellipsoid
    A = (pi / 4 * (b + c) *
         ((b + c) / 2 + 2 * h ** 2 / sqrt(4 * h ** 2 - (b + c) ** 2) * asin(sqrt(4 * h ** 2 - (b + c) ** 2) / (2 * h))))
    if np.isnan(A):
        Warning('sqrt results in negative value since: h<b+c' )
    V = pi / 6 * b * c * h
    D1 = b
    D2 = c
    D3 = h
    return V, A, np.sort([D1, D2, D3])


def shape_5(d, h):  # Cone
    A = pi / 4 * d * (d + sqrt(4 * h ** 2 + d ** 2))
    V = pi / 12 * d ** 2 * h
    D1 = d
    D2 = d
    D3 = h
    return V, A, np.sort([D1, D2, D3])


def shape_7(a, b, c):  # Parallelepiped
    A = 2 * a * b + 2 * b * c + 2 * a * c
    V = a * b * c
    D1 = a
    D2 = b
    D3 = c
    return V, A, np.sort([D1, D2, D3])


def shape_8(a, b, c):  # Prism on elliptic base
    # A   = pi/2 * (a*b + c * (a + b))
    # this formula from http://phytobioimaging.unisalento.it/en-us/products/AtlasOfShapes.aspx?ID_Tipo=3
    # suggests a first order approximation for the Area, better this
    A = c * pi * ((1 / 2) * a + (1 / 2) * b) * (1 + (a - b) ** 2 / (4 * (a + b) ** 2)) + (1 / 2) * pi * a * b
    V = (1 / 4) * pi * a * b * c
    D1 = a
    D2 = b
    D3 = c
    return V, A, np.sort([D1, D2, D3])


def shape_9(a, b, c):  # Prism on parallelogram base
    # A   = a*b + sqrt(a**2 + b**2)/4 * c #this
    # formula from Hillebrand 1999 is wrong, it gives approximately 50# error
    A = a * b + 2 * sqrt(a ** 2 + b ** 2) * c
    V = (1 / 2) * a * b * c
    D1 = a
    D2 = b
    D3 = c
    return V, A, np.sort([D1, D2, D3])


def shape_10(a):  # Cube
    A = 6 * a ** 2
    V = a ** 3
    D1 = a
    D2 = a
    D3 = a
    return V, A, np.sort([D1, D2, D3])


def shape_11(a, c):  # Prism on triangle base 1
    # there is a typo in the formula for area 3 ab + \sqrt(3)/2 *
    # a^2 on web, it should be
    A = 3 * a * c + (1 / 2) * sqrt(3) * a ** 2
    V = (1 / 4) * sqrt(3) * a ** 2 * c
    D1 = a
    D2 = sqrt(3) / 2 * a  # height
    D3 = c
    return V, A, np.sort([D1, D2, D3])


def shape_12(a, b, c):  # Half prism on elliptic base
    # A = pi/4 * (a * b + b * c  + a*c ) + a*c
    # this formula from http://phytobioimaging.unisalento.it/en-us/products/AtlasOfShapes.aspx?ID_Tipo=3
    # is wrong (also in Hillebrand 1999)
    # one can use A = pi/2 * a*b + pi/2 * (a/2 + b)*c + a*c
    # a more precise formula is
    A = (1 / 2) * pi * ((1 / 2) * a + b) * (1 + (a - 2 * b) ** 2 / (4 * (a + 2 * b) ** 2)) * c + a * c + (
                1 / 2) * pi * a * b
    V = (1 / 4) * pi * a * b * c
    D1 = a
    D2 = b
    D3 = c
    return V, A, np.sort([D1, D2, D3])


def shape_14(d, h):  # Double cone
    A = (1 / 2) * pi * d * sqrt(d ** 2 + h ** 2)
    V = (1 / 12) * pi * d ** 2 * h
    D1 = d
    D2 = d
    D3 = h
    return V, A, np.sort([D1, D2, D3])


def shape_15(d1, d2, h):  # 2 truncated cones. I did not check it (there is only one cell)
    l = sqrt((d2 / 2 - d1 / 2) ** 2 + h ** 2)
    A = pi * l * (d1 + d2) + pi / 2 * d1 ** 2
    V = pi / 6 * h * (d1 ** 2 + d1 * d2 + d2 ** 2)
    D1 = max(d1, d2)
    D2 = max(d1, d2)
    D3 = 2 * h
    return V, A, np.sort([D1, D2, D3])


def shape_16(d, d1, d2, h, h1, h2):  # Prolate spheroid + 2 Cylinders
    A = pi * d1 * h1 + pi * d2 * h2 + pi * d * h ** 2 * asin(sqrt(-d ** 2 + h ** 2) / h) / (
                2 * sqrt(-d ** 2 + h ** 2)) + (1 / 2) * pi * d ** 2
    V = (1 / 4) * pi * d1 ** 2 * h1 + (1 / 4) * pi * d2 ** 2 * h2 + (1 / 6) * pi * d ** 2 * h
    D1 = d
    D2 = d
    D3 = h1 + h2 + h
    return V, A, np.sort([D1, D2, D3])


def shape_17(d, h):  # Cylinder +2 cones # I assume that the cones are equilateral
    # Both Hillebrand 1999 and web have mistakes in the formulas
    # A   = pi*d**2+pi*d*(h-sqrt(3)*d)
    # V = (1/12)*pi*d**3*sqrt(3)+(1/4)*pi*d**2*(h-sqrt(3)*d)
    # here I assume that the cone height equals d/2
    A = (1 / 2) * pi * d ** 2 * sqrt(2) + pi * d * (h - d)
    V = -(1 / 6) * pi * d ** 3 + (1 / 4) * pi * d ** 2 * h
    D1 = d
    D2 = d
    D3 = h
    return V, A, np.sort([D1, D2, D3])


def shape_19(d, h):  # cone+half sphere. The Volume is wrong on web, Area is wrong on web
    A = (1 / 4) * pi * d * (sqrt(2 * d ** 2 - 4 * d * h + 4 * h ** 2) + 2 * d)
    V = (1 / 24) * pi * d ** 2 * (d + 2 * h)
    D1 = d
    D2 = d
    D3 = h
    return V, A, np.sort([D1, D2, D3])


def shape_20(b, c, h, h1):  # Half ellipsoid + Cone (on elliptic base)
    # on web there was only Volume
    ConeSideArea = (1 / 2) * pi * ((1 / 4) * b * sqrt(b ** 2 + 4 * h1 ** 2) + (1 / 4) * c * sqrt(c ** 2 + 4 * h1 ** 2))
    HalfEllArea = ((1 / 8) * pi * (b + c) *
                   ((1 / 2) * b + (1 / 2) * c + 2 * h ** 2 * asin(sqrt(4 * h ** 2 - (b + c) ** 2) / (2 * h)) /
                    sqrt(4 * h ** 2 - (b + c) ** 2)))
    A = ConeSideArea + HalfEllArea
    V = (1 / 12) * pi * c * b * h1 + (1 / 12) * pi * b * c * h
    D1 = b
    D2 = c
    D3 = h / 2 + h1
    return V, A, np.sort([D1, D2, D3])


def shape_21(a, a1, b, b1, c):  # Prism on elliptic base+ box, ASh gives a correct formula for V, but a bit wrong for A
    P = pi * ((1 / 2) * a + (1 / 2) * b) * (1 + (a - b) ** 2 / (4 * (a + b) ** 2)) + 2 * a1
    A = P * c + 2 * a1 * b1 + (1 / 2) * pi * a * b
    V = a1 * b1 * c + (1 / 4) * pi * a * b * c
    D1 = b
    D2 = c
    D3 = a + a1
    return V, A, np.sort([D1, D2, D3])


def shape_22(d, h):  # Cylinder + 2 Half spheres
    A = d * pi * (d + h)
    V = (1 / 6) * pi * d ** 3 + (1 / 4) * pi * d ** 2 * h
    D1 = d
    D2 = d
    D3 = h + d
    return V, A, np.sort([D1, D2, D3])


def shape_23(b, c, d1, d2, d3, h, h1, h2, h3):  # Ellipsoid+2cones+cylinder
    A = ((1 / 4) * pi * (b + c) *
         ((1 / 2) * b + (1 / 2) * c + 2 * h ** 2 * asin(sqrt(4 * h ** 2 - (b + c) ** 2) / (2 * h)) / sqrt(
             4 * h ** 2 - (b + c) ** 2))
         - (1 / 4) * pi * d2 ** 2 - (1 / 4) * pi * d3 ** 2 + h1 * d1 * pi + (1 / 2) * pi * d2 * sqrt(
                h2 ** 2 + (1 / 4) * d2 ** 2)
         + (1 / 2) * pi * d3 * sqrt(h3 ** 2 + (1 / 4) * d3 ** 2))
    V = (1 / 6) * pi * b * c * h + (1 / 4) * pi * d1 ** 2 * h1 + (1 / 12) * pi * d2 ** 2 * h2 + (
                1 / 12) * pi * d3 ** 2 * h3
    D1 = b
    D2 = c
    D3 = h + h1 + max(h2, h3)
    return V, A, np.sort([D1, D2, D3])


def shape_24(b, c, d1, h, h1):  # Ellipsoid + Cone
    A = ((1 / 4) * pi * (b + c) *
         ((1 / 2) * b + (1 / 2) * c + 2 * h ** 2 * asin(sqrt(4 * h ** 2 - (b + c) ** 2) / (2 * h)) /
          sqrt(4 * h ** 2 - (b + c) ** 2))
         - (1 / 4) * pi * d1 ** 2 + (1 / 2) * pi * d1 * sqrt(h1 ** 2 + (1 / 4) * d1 ** 2))
    V = (1 / 6) * pi * b * c * h + (1 / 12) * pi * d1 ** 2 * h1
    D1 = b
    D2 = c
    D3 = h + h1
    return V, A, np.sort([D1, D2, D3])


def shape_25(d1, d2, d3, d4, h1, h2, h3, h4):  # Cylinder + 3 Cones
    A = (pi * ((1 / 2) * d1 + (1 / 2) * d4) * sqrt(((1 / 2) * d4 - (1 / 2) * d1) ** 2 + h4 ** 2)
         + (1 / 4) * pi * d4 ** 2 + h1 * d1 * pi + (1 / 4) * pi * d1 ** 2 - (1 / 4) * pi * d3 ** 2 - (
                     1 / 4) * pi * d2 ** 2
         + (1 / 2) * pi * d3 * sqrt(h3 ** 2 + (1 / 4) * d3 ** 2) + (1 / 2) * pi * d2 * sqrt(
                h2 ** 2 + (1 / 4) * d2 ** 2))
    V = ((1 / 12) * pi * h4 * (d1 ** 2 + d1 * d4 + d4 ** 2)
         + (1 / 4) * pi * d1 ** 2 * h1 + (1 / 12) * pi * d2 ** 2 * h2 + (1 / 12) * pi * d3 ** 2 * h3)
    D1 = d4
    D2 = d4
    D3 = h1 + h4 + max(h2, h3)
    return V, A, np.sort([D1, D2, D3])


def shape_27(d):  # Half sphere
    A = 3 * pi * d ** 2 * (1 / 4)
    V = (1 / 12) * pi * d ** 3
    D1 = d
    D2 = d
    D3 = d / 2
    return V, A, np.sort([D1, D2, D3])


def shape_34(a, b1, b2, c, c1, c2, h1, h2):  # 2 half ellipsoids + prism on elliptic base
    A = ((1 / 8) * pi * (2 * b1 + c1) *
         (b1 + (1 / 2) * c1 + 2 * h1 ** 2 * asin(sqrt(4 * h1 ** 2 - (2 * b1 + c1) ** 2) / (2 * h1)) / sqrt(
             4 * h1 ** 2 - (2 * b1 + c1) ** 2))
         + (1 / 4) * pi * h1 * c1 - (1 / 4) * pi * a * c1
         + (1 / 8) * pi * (2 * b2 + c2) *
         (b2 + (1 / 2) * c2 + 2 * h2 ** 2 * asin(sqrt(4 * h2 ** 2 - (2 * b2 + c2) ** 2) / (2 * h2)) / sqrt(
             4 * h2 ** 2 - (2 * b2 + c2) ** 2))
         + (1 / 4) * pi * h2 * c2 - (1 / 4) * pi * a * c2 + pi * ((1 / 2) * a + (1 / 2) * c1) * (
                     1 + (a - c1) ** 2 / (4 * (a + c1) ** 2)) * c)
    # V   = (1/12)*pi*b1*c1*h1+(1/12)*pi*b2*c2*h2+(1/4)*pi*a*b*c
    V = (1 / 6) * pi * b1 * c1 * h1 + (1 / 6) * pi * b2 * c2 * h2 + (1 / 4) * pi * a * c1 * c
    D1 = np.mean([c1, c2], 2)
    D2 = b1 + c + b2
    D3 = np.mean([h1, h2], 2)
    return V, A, np.sort([D1, D2, D3])


def shape_35(a, b, c):  # Cymbelloid. in Hi99 we do not have area, in ASh the area is wrongly found
    # (it should be b instead of c and arcsin(beta))
    A = (b * (2 * b + 2 * a ** 2 * asin(sqrt(4 * a ** 2 - 16 * b ** 2) / (2 * a)) / sqrt(4 * a ** 2 - 16 * b ** 2))
         * asin(c / (2 * b)) + (1 / 2) * pi * a * b)
    V = 2 * a * b ** 2 * asin(c / (2 * b)) * (1 / 3)
    D1 = b
    D2 = c
    D3 = a
    return V, A, np.sort([D1, D2, D3])


def shape_38(d, h):  # Half Cone (WARNING: Not provided in cellgeom)
    # V=pi/6*D^2*G  (Hillebrand99)
    # A=pi*D*l (Hillebrand99)
    l = np.sqrt((h ** 2 + (d / 2) ** 2))
    # l: length of the diagonal connecting the tip of the cone to any point around the circular base
    V = pi / 12 * d ** 2 * h /2
    A = pi / 2 * d * (d / 2 + l) / 2 + d*h/2  # Last term is the triangle
    D1 = d/2
    D2 = l
    D3 = h
    return V, A, np.sort([D1, D2, D3])


def shape_40(b, c, h):  # Gomphonemoid
    A = (1 / 2) * b * (2 * h + pi * h * asin(c / (2 * h)) + ((1 / 2) * pi - 2) * b)
    V = (1 / 4) * h * b * (h + ((1 / 4) * pi - 1) * b) * asin(c / (2 * h))
    D1 = b
    D2 = c
    D3 = h
    return V, A, np.sort([D1, D2, D3])


def shape_41(a, c, h):  # Sickle-shaped prism
    b = a
    b2 = 0 * b  # assume the inner semi axis equals 0
    A = (1 / 2) * pi * (b * c + b * h + b2 * c - b2 * h + c * h)
    V = (1 / 4) * pi * c * h * a
    D1 = a
    D2 = c
    D3 = h
    return V, A, np.sort([D1, D2, D3])


def shape_43(a, b, c, d1, d2, d3, d4, h1, h2, h3, h4):  # Prism on elliptic base + 4 Cones
    A = (c * pi * ((1. / 2.) * a + (1. / 2.) * b) * (1 + (a - b) ** 2. / (4 * (a + b) ** 2))
         + (1. / 2.) * pi * a * b - (1. / 4.) * pi * d1 ** 2 - (1 / 4) * pi * d2 ** 2 - (1 / 4) * pi * d3 ** 2 - (
                     1 / 4) * pi * d4 ** 2
         + (1 / 2) * pi * d1 * sqrt(h1 ** 2 + (1 / 4) * d1 ** 2) + (1 / 2) * pi * d2 * sqrt(h2 ** 2 + (1 / 4) * d2 ** 2)
         + (1 / 2) * pi * d3 * sqrt(h3 ** 2 + (1 / 4) * d3 ** 2) + (1 / 2) * pi * d4 * sqrt(
                h4 ** 2 + (1 / 4) * d4 ** 2))
    V = ((1 / 4) * pi * a * b * c + (1 / 12) * pi * d1 ** 2 * h1 + (1 / 12) * pi * d2 ** 2 * h2
         + (1 / 12) * pi * d3 ** 2 * h3 + (1 / 12) * pi * d4 ** 2 * h4)
    D1 = a
    D2 = b
    D3 = max(h1, h2) + c + max(h3, h4)
    return V, A, np.sort([D1, D2, D3])


def shape_44(d, h):  # Pyramid (rectangular base)
    A = sqrt(d ** 2 + 4 * h ** 2) * d + d ** 2
    V = (1 / 3) * d ** 2 * h
    D1 = d
    D2 = d
    D3 = h
    return V, A, np.sort([D1, D2, D3])


def shape_46(a, b):  # Prism on triangle-base 2
    A = 3 * a * b + (1 / 2) * a ** 2 * sqrt(3)
    V = (1 / 4) * a ** 2 * sqrt(3) * b
    D1 = a
    D2 = sqrt(3) / 2 * a
    D3 = b
    return V, A, np.sort([D1, D2, D3])


def shape_51(b1, b2, c1, c2, h1, h2):  # 2 Half ellipsoids
    A = ((1 / 8) * pi * (2 * b1 + c1) *
         (b1 + (1 / 2) * c1 + 2 * h1 ** 2 * asin(sqrt(4 * h1 ** 2 - (2 * b1 + c1) ** 2) / (2 * h1)) / sqrt(
             4 * h1 ** 2 - (2 * b1 + c1) ** 2))
         + (1 / 8) * pi * (2 * b2 + c2) *
         (b2 + (1 / 2) * c2 + 2 * h2 ** 2 * asin(sqrt(4 * h2 ** 2 - (2 * b2 + c2) ** 2) / (2 * h2)) / sqrt(
             4 * h2 ** 2 - (2 * b2 + c2) ** 2)))
    V = (1 / 6) * pi * b1 * c1 * h1 + (1 / 6) * pi * b2 * c2 * h2
    D1 = np.mean([c1, c2], 2)
    D2 = b1 + b2
    D3 = np.mean([h1, h2], 2)
    return V, A, np.sort([D1, D2, D3])


def get_dim(fields,headersin,dim,dataset='Olenina'):
    col = -1
    #try to find the dimension
    if dataset=='Olenina':
        if dim=='L1':
            for coli,colname in enumerate(headersin):
                if colname[0:10]=='Length(l1)':
                    col=coli
        if dim=='L2':
            for coli,colname in enumerate(headersin):
                if colname[0:10]=='Length(l2)':
                    col=coli
        if dim=='W':
            for coli,colname in enumerate(headersin):
                if colname[0:8]=='Width(w)':
                    col=coli
        if dim=='H':
            for coli,colname in enumerate(headersin):
                if colname[0:9]=='Height(h)':
                    col=coli
        if dim=='D1':
            for coli,colname in enumerate(headersin):
                if colname[0:12]=='Diameter(d1)':
                    col=coli
        if dim=='D2':
            for coli,colname in enumerate(headersin):
                if colname[0:12]=='Diameter(d2)':
                    col=coli
    elif dataset == 'Roselli':
        raise(ValueError('Extracting dimensions from Roselli datasets not yet implemented'))

    #if found, extract the string and attempt converting to float
    if col == -1: #
        print ('dimension not found: %s' % (dim))
        val = np.nan
    else:
        valstr=fields[col]
    try:
        val=float(valstr)
    except:
        val=np.nan
    return val