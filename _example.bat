rem first step: create in an empty folder a new image

upperbound-rc anew lenm=12 lent=11

rem compute refinements

upperbound-rc lenm=12 lent=11 p123
upperbound-rc lenm=13 lent=11 p123 save=raw

rem if temporarily renaming *.DATA to *.BMP, files can be viewed

pause
