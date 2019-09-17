# 2019-CodeContest-DetectingCorners
Given image find out  all the corners, which match line length and angle of intersections.


PROBLEM STATEMENT:

Detect corners:

 Given a binary image (https://en.wikipedia.org/wiki/Binary_image) findout
 all the corners.A corner is defined as the intersection of two lines.
 Please use the attached example images.

 Input parameters: <binary_image_path> <minimum_line_length> <angle_min> <angle_max>

 binary_image_path   - Input image on which corner detection needs to be performed.
 minimum_line_length - Length of line segment in pixels
 angle_min           - Minimum angle in degrees
 angle_max           - Maximum angle in degrees

 Once all corners are found:

 1. Print on screen or mark them on the result image.
 2. Print the length of line segments and the angle between the line segments.
 3. If the detected line segment is less than the specified minimum length,
   ignore the corner.
 4. Consider the corner only if the angle between the line segments is
   in between the specified angle range

 5. In the given binary image the line may not be single pixel wide. You
   would need to do pre-processing for getting a single pixel wide line.

 6. Use of opencv or any other image processing libraries are permitted.
   You could build / use any machine learning models as well.

 Success criteria:

 1. Correctness: Program should find all the corners which meets the
   criteria (length of line segment and the angle bweteen the line segment
   is within the specified range)

 2. Speed

 3. Memory usage
