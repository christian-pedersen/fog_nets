import sys


"""
    Script takes radius[0], x_dist[1] and y_dist[2] as commmand line arguments
    and writes a .geo file with 5 circles in a rectangluar domain
    with 2 circles in the front row and 1 + 2*1/2 circles in the back row
"""

arg = sys.argv[1::]

radius = float(arg[0])
x_dist = float(arg[1])
y_dist = float(arg[2])

scale = min([radius, x_dist, y_dist])

lm = scale / 7.5
lc = scale / 15.

if x_dist < 0.1*radius:
    x_dist = 0.1*radius
if y_dist < 0.1*radius:
    y_dist = 0.1*radius

geo_file = open('mesh.geo', 'w')

geo_file.write('SetFactory("OpenCASCADE");\n')

geo_file.write('Rectangle(1) = {0, 0, 0, %f, %f};\n' % ((40*radius + x_dist), (4*radius+2*y_dist)))

xcenter, ycenter = [], []

for i in range(5):
    if i < 2:
        x_center = 5*radius
        y_center = (2*i+1)*radius + 0.5*(2*i+1)*y_dist

    else:
        x_center = (5+2)*radius + x_dist
        y_center = 2*(i-2)*radius + (i-2)*y_dist

    xcenter.append(x_center)
    ycenter.append(y_center)

    geo_file.write('Circle(%g) = {%f, %f, 0, %f, 0, 2*Pi}; Line Loop(%g) = {%g}; Plane Surface(%g) = {%g};\n' % (i+5, x_center, y_center, radius, i+5, i+5, i+5, i+5))

xcenter[2], xcenter[3] = xcenter[3], xcenter[2]
ycenter[2], ycenter[3] = ycenter[3], ycenter[2]

geo_file.write('BooleanDifference(10) = { Surface{1}; Delete; }{ Surface{5, 6, 7, 8, 9}; Delete; };\n')


point_no = 30
cylinder_points = [1, 2, 5, 6, 9, 10, 11]
for i in range(5):
    if i < 3:
        geo_file.write('Point(%g) = {%f, %f, 0}; Point(%g) = {%f, %f, 0}; Point(%g) = {%f, %f, 0};\n' % (point_no+3*i, xcenter[i]-radius, ycenter[i], point_no+(3*i +1), xcenter[i], ycenter[i]+radius, point_no+(3*i +2), xcenter[i], ycenter[i]-radius))
        for j in range(3):
            cylinder_points.append(point_no+3*i + j)
    elif i == 3:
        geo_file.write('Point(%g) = {%f, %f, 0};\n' % (point_no+3*i, xcenter[i], ycenter[i]+radius))
        point_no += 3*i
        cylinder_points.append(point_no)
    else:
        geo_file.write('Point(%g) = {%f, %f, 0};\n' % (point_no+1, xcenter[i], ycenter[i]-radius))
        cylinder_points.append(point_no+1)

geo_file.write('Characteristic Length {3, 4, 7, 8} = %f;\n' % lm)
geo_file.write('Characteristic Length {%g' % cylinder_points[0])

for i in range(len(cylinder_points)-1):
    geo_file.write(', %g' % cylinder_points[i+1])
geo_file.write('} = %f;\n'%lc)

line_points = [3, 2, 5, 4, 1, 8, 7, 6]
surface  =  [4, 2, 6, 8]
surf_no = 50
for i in range(2):
    #geo_file.write('Line(%g) = {%g, %g}; Line(%g) = {%g, %g};\n' % (surf_no+2*i, line_points[4*i], line_points[4*i+1], surf_no+2*i+1, line_points[4*i+2], line_points[4*i+3]))
    #geo_file.write('Plane Surface(%g) = {%g}; Plane Surface(%g) = {%g};\n' % (surf_no+2*i, surface[2*i], surf_no+2*i+1, surface[2*i+1]))
    geo_file.write('Periodic Line {%g} = {%g} Translate {0, %f, 0};\n' % (surface[2*i], surface[2*i+1], 4*radius+2*y_dist))


geo_file.close()
