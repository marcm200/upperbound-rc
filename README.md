# upperbound-rc

Reliable escape time method to describe the exterior of the quadratic Mandelbrot set
using interval and affine arithmetics and outward rounding. 

The orbit for a complex interval is followed by interval and subsequently by affine
arithmetics.

see https://fractalforums.org/index.php?topic=3467.msg27794#msg27794

Parameters: see short helppage calling the software without command line parameters.

# Example usage:

upperbound-rc lenm=12 lent=11 anew load=raw save=raw
upperbound-rc lenm=12 lent=11 p123 load=raw save=raw
upperbound-rc lenm=13 lent=11 p123 load=raw save=raw


# Used flow for the results in subfolder:

upperbound-rc lenm=16 lent=15 anew
upperbound-rc lenm=16 lent=15 p123
upperbound-rc lenm=17 lent=15 p123
upperbound-rc lenm=18 lent=15
upperbound-rc lenm=18 lent=15 skipless=2500 
upperbound-rc lenm=20 lent=15 save=raw skipless=2500




