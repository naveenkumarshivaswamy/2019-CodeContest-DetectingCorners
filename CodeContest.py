import cv2
import intersectModel_point
import math
import numpy as np
import os

# import imutils
# from __future__ import division
# from __future__ import print_function
# import warnings
# warnings.simplefilter(action = "ignore", category = RuntimeWarning)

modelOutput = './ModelOutput/'

def getPointXY(point):
    x, y = point
    # print ("x: ", x)
    # print("y: ", y)
    return x, y

def getLineLength(LineSegment):
    # print("LineSegment: ", LineSegment)
    # print("PointA: ", LineSegment[0])
    x1, y1 = getPointXY(LineSegment[0])
    # print("x1: ", x1)
    # print("y1: ", y1)
    # print("PointB: ", LineSegment[1])
    x2, y2 = getPointXY(LineSegment[1])
    # print("x2: ", x2)
    # print("y2: ", y2)
    
    length = math.sqrt(((x2 - x1) ** 2) + ((y2 - y1) ** 2))
    # # length = math.hypot(x2 - x1, y2 - y1)
    # print("Length of LineSegment: ", length)
    return length

def slope(x1, y1, x2, y2):
    return (y2-y1)/(x2-x1)

def angle(s1, s2):
    return math.degrees(math.atan((s2-s1)/(1+(s2*s1))))

def getAngle(lineA, lineB):
	slope1 = slope(lineA[0][0], lineA[0][1], lineA[1][0], lineA[1][1])
	slope2 = slope(lineB[0][0], lineB[0][1], lineB[1][0], lineB[1][1])
	ang = angle(slope1, slope2)
	# print('Angle in degrees = ', ang)
	return ang

def getStatistics(intersections):
    for interSectionPoint in intersections:
        print ("InterSectionPoint raw data: ", interSectionPoint[0], " --> ", interSectionPoint[1])

        print ("InterSectionPoint: ", interSectionPoint[0])
        x, y = getPointXY(interSectionPoint[0])
        print ("x: ", x)
        print("y: ", y)

        print("LineSegment-1: ", interSectionPoint[1][0])
        print("LineSegment-1 length: ", getLineLength(interSectionPoint[1][0]))
        x1, y1 = getPointXY(interSectionPoint[1][0][0])
        print ("x1: ", x1)
        print("y1: ", y1)
        x2, y2 = getPointXY(interSectionPoint[1][0][1])
        print("x2: ", x2)
        print("y2: ", y2)

        print("LineSegment-2: ", interSectionPoint[1][1])
        print("LineSegment-2 length: ", getLineLength(interSectionPoint[1][1]))
        x1, y1 = getPointXY(interSectionPoint[1][1][0])
        print("x1: ", x1)
        print("y1: ", y1)
        x2, y2 = getPointXY(interSectionPoint[1][1][1])
        print("x2: ", x2)
        print("y2: ", y2)

        print('Angle in degrees = ', getAngle(interSectionPoint[1][0], interSectionPoint[1][1]))
# ----------------------------------------------------------------------------------------------------------------------
# lines will be vertical and horizontal,
# We can easily split the lines based on their position.
# If the two y-coordinates of one line are near to each other, then
#    the line is mostly horizontal.
# If the two x-coordinates of one line are near to each other, then
#    the line is mostly vertical.
# now we segment your lines into vertical lines and horizontal lines
def segmentTheline(imgRaw, lines, delta):
    h_lines = []
    v_lines = []
    for line in lines:
        for x1, y1, x2, y2 in line:
            # Vertical line because x-values are near to delta value
            if abs(x2-x1) < delta:
                v_lines.append(line)
            # Horizontal line because y-values are near to delta value
            if abs(y2-y1) < delta:
                h_lines.append(line)

    # Draw Segmented-Hough-Lines
    houghimg = imgRaw
    for line in h_lines:
        for x1, y1, x2, y2 in line:
            color = [0, 0, 255]  # color hoz lines red
            cv2.line(houghimg, (x1, y1), (x2, y2), color=color, thickness=1)
    for line in v_lines:
        for x1, y1, x2, y2 in line:
            color = [255, 0, 0]  # color vert lines blue
            cv2.line(houghimg, (x1, y1), (x2, y2), color=color, thickness=1)
            
    return h_lines, v_lines, houghimg

def maskBlackColorFromInputImage(imgRaw, fileName):
    # find all the 'black' shapes in the image
    lower = np.array([0, 0, 0])
    upper = np.array([15, 15, 15])
    blackMaskedImg = cv2.inRange(imgRaw, lower, upper)
    # cv2.imshow(fileName+"-maskBlackColorFromInputImage", blackMaskedImg)
    # cv2.waitKey(0)
    # cv2.imwrite(modelOutput + fileName + '-maskBlackColorFromInputImage.png', blackMaskedImg)
    return blackMaskedImg

def getEdgesFromCanny(imgRaw, low_threshold, high_threshold):
    edges = cv2.Canny(imgRaw, low_threshold, high_threshold)
    # edges = cv2.Canny(imgRaw, low_threshold, low_threshold, apertureSize=5)
    # edges = cv2.Canny(imgRaw, low_threshold, low_threshold, None, 3)
    
    # vis = imgRaw.copy()
    # vis = np.uint8(vis/2.)
    # vis[edges != 0] = (0, 255, 0)
    
    # cv2.imshow('edge', vis)
    # cv2.waitKey(0)
    # cv2.imwrite('./a.png', vis)
    return edges
    
def main(binary_image_path, minimum_line_length, angle_min, angle_max):
    fileName = os.path.basename(binary_image_path)
    print ("DEBUG: Processing input image --> " + fileName + " minimum_line_length= " + str(minimum_line_length) + " angle_min= " + str(angle_min) + " angle_max= " + str(angle_max))
    # ----------------------------------------------------------------------------------------------------------------------
    img = cv2.imread(binary_image_path)
    # img = cv2.resize(img, None, fx=.5, fy=.5)
    # img = maskBlackColorFromInputImage(img.copy(), fileName)
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    blur_gray = cv2.GaussianBlur(gray, (5, 5), 0)
    edges = getEdgesFromCanny(blur_gray, 50, 150)
    # ----------------------------------------------------------------------------------------------------------------------
    # Running the probabilistic HoughLines transformation [to get 4 lines around the square so that I can obtain the vertices of the square]
    # Making a grid of lines by passing through ton of angles and distances.
    # If lines go over any white pixels from Canny-dilated-image, then it gives that line a SCORE for each point it goes through.
    # output contain line segments, where all lines will have 2 points.
    # Run Hough on edge detected image
    # Output "lines" is an array containing endpoints of detected line segments
    lines = []
    # lines = cv2.HoughLines(edges, 1, np.pi / 180, 15, np.array([]), min_theta = angle_min, max_theta = angle_max)
    lines = cv2.HoughLinesP(edges, 1, np.pi / 180 , 15, np.array([]), minimum_line_length, 20)
    # print("DEBUG: HoughLines--> ", lines)
    
    points = []
    line_image = np.copy(img) * 0
    for line in lines:
        for x1, y1, x2, y2 in line:
            points.append(((x1 + 0.0, y1 + 0.0), (x2 + 0.0, y2 + 0.0)))
            cv2.line(line_image, (x1, y1), (x2, y2), (0, 255, 0), 2)
    lines_edges = cv2.addWeighted(img, 0.8, line_image, 1, 0)
    # print("DEBUG: lines_edges --> ", lines_edges.shape)
    
    # cv2.imshow(fileName+"-Hough-Lines-image", lines_edges)
    # cv2.waitKey(0)
    # cv2.imwrite(modelOutput+ fileName+'-Hough-Lines-image.png', lines_edges)
    # -----------------------------------------------------------------------------------------------------------------
    # segment the lines to 2 classes, one horizontal & one vertical
    delta = 150
    h_lines, v_lines, houghimg = segmentTheline(img,lines, delta)
    
    # cv2.imshow(fileName+"-Segmented-Hough-Lines-image", houghimg)
    # cv2.waitKey(0)
    # cv2.imwrite(modelOutput + fileName + '-Segmented-Hough-Lines-image.png', houghimg)
    # -----------------------------------------------------------------------------------------------------------------
    # print ("DEBUG: Points--> ", points)
    
    # # noteA --> Only intersection points
    # intersections = intersectModel_point.isect_segments(points)
    # print ("DEBUG: Intersection points--> ", intersections)
    # # Strategy to remove multiple intersection points in a small area:
    # for idx, inter in enumerate(intersections):
    #     a, b = inter
    #     match = 0
    #     for other_inter in intersections[idx:]:
    #         c, d = other_inter
    #         if abs(c-a) < 15 and abs(d-b) < 15:
    #             match = 1
    #             intersections[idx] = ((c+a)/2, (d+b)/2)
    #             intersections.remove(other_inter)
    #     if match == 0:
    #         intersections.remove(inter)
    # for interSectionPoint in intersections:
    #     a, b = interSectionPoint
    #     for i in range(10):
    #         for j in range(10):
    #             lines_edges[int(b) + i, int(a) + j] = [0, 0, 255]

    # noteB --> Intersection points with corresponding line segments
    intersections = intersectModel_point.isect_segments_include_segments(points)
    # print ("DEBUG: Intersection points--> ", intersections)
    print ("DEBUG: Intersection points--> ", getStatistics(intersections))
    lines_intersection = lines_edges.copy()
    for interSectionPoint in intersections:
        a, b = interSectionPoint[0]
        for i in range(10):
            for j in range(10):
                lines_intersection[int(b) + i, int(a) + j] = [0, 0, 255]
    # ----------------------------------------------------------------------------------------------------------------------
    # cv2.imshow(fileName+"-Intersection-Points-image", lines_intersection)
    # cv2.waitKey(0)
    # cv2.imwrite(modelOutput+ fileName+'-Intersection-Points-image.png', lines_intersection)
    
    #To view all images in oneshot axis=1=horizantal axis=0=verticle
    img_concatenate = np.concatenate((lines_edges,houghimg,lines_intersection),axis=1)
    # cv2.imshow('SeriesOfImageProcessing', img_concatenate)
    cv2.imwrite(modelOutput + fileName + '-Intersection-Points-image.png', img_concatenate)
    # cv2.waitKey(0)
    # # ----------------------------------------------------------------------------------------------------------------------
    
if __name__ == '__main__':
    from glob import glob
    
    # for binary_image_path in glob('./TestDataSet/Test*'):
    for binary_image_path in glob('./ProblemStatement/edge*.png'):
        # Input parameters: <binary_image_path> <minimum_line_length> <angle_min> <angle_max>
        main(binary_image_path, minimum_line_length = 50, angle_min = 1, angle_max = 180)
    cv2.destroyAllWindows()
    print('Done')
    
    # # ----------------------------------------------------------------------------------------------------------------
    # # Reading input manually
    # # Input parameters: <binary_image_path> <minimum_line_length> <angle_min> <angle_max>
    # binary_image_path = './TestDataSet/Test1.png'
    # binary_image_path = './TestDataSet/Test2.jpg'
    # binary_image_path = './TestDataSet/Test3.png'
    # binary_image_path = './TestDataSet/Test4.jpg'
    # binary_image_path = './TestDataSet/Test5.jpg'
    # binary_image_path = './TestDataSet/Test5.png'
    # binary_image_path = './TestDataSet/Test7.png'
    # binary_image_path = './TestDataSet/Test8.png'
    #
    # binary_image_path = './ProblemStatement/edge.png'
    # binary_image_path = './ProblemStatement/edge_2.png'
    # binary_image_path = './ProblemStatement/edge_3.png'
    # binary_image_path = './ProblemStatement/edge_4.png'
    #
    # main(binary_image_path, minimum_line_length=50, angle_min=0, angle_max=0)
    # cv2.destroyAllWindows()
    # print('Done')
    # # ----------------------------------------------------------------------------------------------------------------