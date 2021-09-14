import numpy as np

""" Preprocessing """
# Read the turbine's operational parameters
dataFile = open('data', 'r')
ri = float( dataFile.readline().rstrip(';\n').split()[1] )
ro = float( dataFile.readline().rstrip(';\n').split()[1] )
T =  float( dataFile.readline().rstrip(';\n').split()[1] )
Q =  float( dataFile.readline().rstrip(';\n').split()[1] )
rho =float( dataFile.readline().rstrip(';\n').split()[1] )
dataFile.close()

# Read the coordinates
points = np.loadtxt('points.csv', delimiter=',')
start = []
end = []
for i in range(points.shape[0]):
    start.append(list(points[i ,0:3]))
    end.append(list(points[i ,3:]))

# determine limits of the domain:
x = np.hstack((points[:,0], points[:,3]))
z = np.hstack((points[:,2], points[:,5]))

xmin = x.min() - 5.00 * ro
xmax = x.max() + 25.0 * ro
zmin = z.min() - 5.00 * ro
zmax = z.max() + 5.00 * ro
ymin = -5.00 * ro
ymax =  5.00 * ro
pom = xmin + 0.5
ny = 20
nz = 20
nx = ny * (xmax - xmin) / (ymax - ymin)

""" Domain and mesh file """
# write the "domain" file:
domainFile = open('domain', 'a')
domainFile.write('xmin \t\t {0:7.3f};\n'.format(xmin))
domainFile.write('xmax \t\t {0:7.3f};\n'.format(xmax))
domainFile.write('zmin \t\t {0:7.3f};\n'.format(zmin))
domainFile.write('zmax \t\t {0:7.3f};\n'.format(zmax))
domainFile.write('ymin \t\t {0:7.3f};\n'.format(ymin))
domainFile.write('ymax \t\t {0:7.3f};\n'.format(ymax))
domainFile.write('POM  \t\t {0:7.3f};\n'.format(pom))
domainFile.write('nx  \t\t {0:7d};\n'.format( int(nx) ))
domainFile.write('ny  \t\t {0:7d};\n'.format( int(ny) ))
domainFile.write('nz  \t\t {0:7d};\n'.format( int(nz) ))
domainFile.close()

""" SnappyHexMesh geometries """
# number of disks:
n = points.shape[0]

template3 = '\tcylinder3{0:}\n\
            \t{{\n\
            \t\ttype   \t\t searchableCylinder;\n\
            \t\tpoint1 \t\t ({1:7.3f} {2:7.3f} {3:7.3f});\n\
            \t\tpoint2 \t\t ({4:7.3f} {5:7.3f} {6:7.3f});\n\
            \t\tradius \t\t {7:7.3f};\n\
            \t}}'

template2 = '\tcylinder2{0:}\n\
            \t{{\n\
            \t\ttype   \t\t searchableCylinder;\n\
            \t\tpoint1 \t\t ({1:7.3f} {2:7.3f} {3:7.3f});\n\
            \t\tpoint2 \t\t ({4:7.3f} {5:7.3f} {6:7.3f});\n\
            \t\tradius \t\t {7:7.3f};\n\
            \t}}'

template1 = '\tcylinder1{0:}\n\
            \t{{\n\
            \t\ttype   \t\t searchableCylinder;\n\
            \t\tpoint1 \t\t ({1:7.3f} {2:7.3f} {3:7.3f});\n\
            \t\tpoint2 \t\t ({4:7.3f} {5:7.3f} {6:7.3f});\n\
            \t\tradius \t\t {7:7.3f};\n\
            \t}}'

geoFile = open('geometry', 'w')

offset3 = 0.15
offset2 = 0.50
offset1 = 1.00
radius3 = ro * (1 + offset3)
radius2 = ro * (1 + offset2)
radius1 = ro * (1 + offset1)

for obj in range(n):
    xstart = end[obj][0] 
    xend = start[obj][0]
    geoFile.write(template3.format(obj,
                              xstart - ro,
                              end[obj][1],
                              end[obj][2],
                              xend + ro,
                              start[obj][1],
                              start[obj][2],
                              radius3))
    geoFile.write('\n\n')
    geoFile.write(template2.format(obj,
                              xstart - (2.0 * ro),
                              end[obj][1],
                              end[obj][2],
                              xend + (4.0 * ro),
                              start[obj][1],
                              start[obj][2],
                              radius2))
    geoFile.write('\n\n')
    geoFile.write(template1.format(obj,
                              xstart - (3.0 * ro),
                              end[obj][1],
                              end[obj][2],
                              xend + (6.0 * ro),
                              start[obj][1],
                              start[obj][2],
                              radius1))
    geoFile.write('\n\n')
    # end of loop
    
geoFile.close()
    
""" SnappyHexMesh refinementRegions """
regionFile = open('regions' ,'w')

region3 = '\tcylinder3{0:}\n\
        \t{{\n\
        \t\tmode \t\t inside;\n\
        \t\tlevels \t\t ((3 3));\n\
        \t}}'

region2 = '\tcylinder2{0:}\n\
        \t{{\n\
        \t\tmode \t\t inside;\n\
        \t\tlevels \t\t ((2 2));\n\
        \t}}'

region1 = '\tcylinder1{0:}\n\
        \t{{\n\
        \t\tmode \t\t inside;\n\
        \t\tlevels \t\t ((1 1));\n\
        \t}}'
        
for obj in range(n):
    regionFile.write(region3.format(obj))
    regionFile.write('\n\n')
    regionFile.write(region2.format(obj))
    regionFile.write('\n\n')
    regionFile.write(region1.format(obj))
    regionFile.write('\n\n')
    # end of loop

regionFile.close()

""" Merge snappyHexMesh files """
filenames = ['bit1', 'geometry', 'bit2', 'regions', 'bit3']
with open('snappyHexMeshDict', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            outfile.write(infile.read())

""" fvSolution """
from math import pi

# Write the fvSolution variable file
variableFile = open('fvSolutionVariablePart', 'w')

variableFile.write('numberOfDisks\n')
variableFile.write('{\n')
variableFile.write('\tnod\t\t{};\n'.format(len(start)))
variableFile.write('}\n\n')

for ix, (sp, ep) in enumerate( zip(start, end) ):
    variableFile.write('actuatorDisk{}\n'.format(ix)) 
    variableFile.write('{\n')
    variableFile.write('\tinteriorRadius\t{0:.2f};\n'.format(ri))
    variableFile.write('\texteriorRadius\t{0:.2f};\n'.format(ro))
    variableFile.write('\tthrust\t\t\t{0:.2f};\n'.format(T))
    variableFile.write('\ttorque\t\t\t{0:.2f};\n'.format(Q))
    variableFile.write('\tdensity\t\t\t{0:.2f};\n'.format(rho))
    variableFile.write('\tstartPoint\t\t({0:.2f} {1:.2f} {2:.2f});\n'.format(*sp))
    variableFile.write('\tendPoint\t\t({0:.2f} {1:.2f} {2:.2f});\n'.format(*ep))
    variableFile.write('}\n\n')

variableFile.close()

# Link the fvSolution files
filenames = ['fvSolutionSettings', 'fvSolutionVariablePart']
with open('fvSolution', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            outfile.write(infile.read())





