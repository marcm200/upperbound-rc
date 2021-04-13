Mandelbrot set z^2+c in a reliable manner (via interval
and affine arithemtics, outward rounding). A pixel
represents a complex intervcal of width 2^-18.

Pixel color denotes the fate of all complex numbers inside a pixel.

White = definite in the exteiror of the Mset.
Gray = unjkudged so far. Contains a superset of the Mset itself

Files to describe an image of size 2^20 x 2^20 pixels depicting the 
complex range [-2..2] + i*[-2..2].

Screen pixel coordinates start at 0,0 at the lower left and
go upwards and rightwards.

Only the lower half is stored.

The whole image was partitioned into chunks of 2^15 x 2^15 pixels.

Three types of files are present:

type 1: fully white

_cert_M0001048576_y0000098304_x0000098304_fullywhite

denotes a chunk that only consists of white (exterior) pixels
and is not stored as an image itself but just as a text file.
y,x denote the screen coordinates of the lower left pixel of the chunk
in global coordinates.


 type 2: fully gray = unjudged

_cert_M0001048576_y0000393216_x0000491520_fullygray


type 3: image both containing exterior and unjudged

_cert_M0001048576_y0000393216_x0000327680.bmp

stored as an 8-bit bitmap with two colors: white and gray.
see e.g. https://fractalforums.org/index.php?topic=3467.msg27794#msg27794


Marc Meidlinger
marc.meidlinger@web.de
April 2021





