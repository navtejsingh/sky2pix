#!/usr/bin/env python

"""
    ------------------------------------------------------------------------
    Routine to translate RA and DEC sky coordinates to X, Y pixel values in
    multiprocessing MPI mode

    Usage: mpiexec -np <ncpus> python sky2pix.py in.fits in_rd_deg.cat --quiet

    Input:
        ncpus: Number of physical cores or processers (or processes)
        image: Input image with basic header keywords
        skyfile: Input file with (ra, dec) pixel value in column 1 and 2
      
      
    [Options]:
        --help: help
        --version: program version
        --verbose: show result messages
        --quiet: don't show result messages
        --filename: output file name (default is pix2sky.dat)
        --degree: (ra, dec) in degrees? (default is yes)
        
        
    Output:
        sky2pix.out: Output file with (ra, dec) and corresponding (x, y)
        pixel values
        
        
    Author:
        Navtej Singh


    Organization:
        Centre for Astronomy, National University of Ireland, Galway, Ireland


    Version:
        10 January 2012     1.0     Original version
    ------------------------------------------------------------------------            
"""

# Load python modules to be used in the routine
from optparse import OptionParser
import sky2pix_mpi as sm
    
    
# Entry point for SKY2PIX utility
if __name__ == '__main__':
    usage = "Usage: mpiexec -np ncpus python %prog [options] image skyfile"
    description = "Description. Utility to convert RA/DEC sky coordinates to X/Y pixel values in multiprocessing MPI mode."
    parser = OptionParser(usage = usage, version = "%prog 1.0", description = description)
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default = False,
                      help = "print result messages to stdout"
                      )
    parser.add_option("-q", "--quiet",
                    action="store_false", dest="verbose", default = True,
                    help = "don't print result messages to stdout"
                    )
    parser.add_option("-d", "--degree", dest = "degree", metavar="DEGREE",
                    action="store", help = "ra/dec in degree? [default is yes]",
                    choices=['yes', 'no'], default = 'yes'
                    )    
    parser.add_option("-f", "--filename", dest = "filename",
                    action='store', metavar="FILE", help = "output file name [default is sky2pix.out]"
                    )
    (options, args) = parser.parse_args()
    
    # Check for number of input arguments
    if len( args ) != 2:
        parser.error("Incorrect number of arguments. Check help for more details.")

    sm.sky2pix(args[0], args[1], options.degree, options.filename, options.verbose)
    
