Varius bits of code playing around with various window functions and some signal processing stuff.
Lots of the code is pretty bad/broken so don't take any of it too seriously it's a bit of a mess, check over everything before using any of this in production.
Some sections had a bit of AI use (mostly the DPSS stuff).

The main interesting things here are:

DPSS window using faer to generate coefficients.

Fractional kaiser window using the fractional modified Bessel function of the first kind mentioned in these papers:

K. Avci, "Fractional Kaiser Window With Application to Finite Impulse Response Digital Filter Design,"
IEEE Access, vol. 12, pp. 155549-155563, 2024
doi: 10.1109/ACCESS.2024.3484930
https://ieeexplore.ieee.org/document/10731688

Martín, A.; Estrada, E. "Fractional-Modified Bessel Function of the First Kind of Integer Order"
MDPI Mathematics 2023, 11, 1630.
doi: 10.3390/math11071630 
https://www.mdpi.com/2227-7390/11/7/1630
