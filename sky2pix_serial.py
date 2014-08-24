#!/usr/bin/env python

"""
    ------------------------------------------------------------------------
    Routine to translate RA and DEC sky coordinates to X, Y pixel values in
    serial mode

    Usage: python sky2pix_serial.py [options] image skyfile

    Input:
        image: Input image with basic header keywords
        skyfile: Input file with (ra, dec) sky value in column 1 and 2
      
      
    [Options]:
        --help: help
        --version: program version
        --verbose: show result messages
        --quiet: don't show result messages
        --filename: output file name (default is sky2pix.out)
        --hour: input (ra, dec) in hms? (default is yes)
        
        
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
from os.path import join, exists
from sys import stderr, stdout, exit
from StringIO import StringIO
from math import cos, sin, radians, degrees
from optparse import OptionParser


# Get header keyword values
def getHeader(image):
    print >> stdout, '\n Getting image header keywords...'
    
    # Input can be a single FITS image, multi-extension image or
    # multi-extension image with particular extension
    if len(image.split('[', 1)) > 1:
        ext = image.split('[', 1)[1].replace(']', '')
        image = image.split('[', 1)[0]
    else:
        ext = ''
    
    # Open the Header Unit List (HDU) to read header keywords
    try:
        hdulist = pyfits.open(image)
    except:
        print >> stderr, 'Error: Not able to read FITS header. Exiting.'
        exit(-1)

    hdulist.close()

    # Get header parameters - checking number of extensions and using 1st extension
    # in case of multi extension FITS image
    if len(hdulist) > 1:
        if ext == '':
            hdrdata = hdulist[1].header
        else:
            
            hdrdata = hdulist[int(ext)].header 
    else:
        hdrdata = hdulist[0].header

    # Get CRPIX keyword values
    crpix1 = hdrdata['CRPIX1']
    crpix2 = hdrdata['CRPIX2']

    # Get CRVAL keyword values
    ra0 = hdrdata['CRVAL1']
    dec0 = hdrdata['CRVAL2']

    # Get CD keyword values
    cd11 = hdrdata['CD1_1']
    cd12 = hdrdata['CD1_2']
    cd21 = hdrdata['CD2_1']
    cd22 = hdrdata['CD2_2']

    # Return keyword values
    return ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22


# Convert RA,DEC from hour,min,sec to degree
def hms2degree(hms):
    return int(hms.split(':')[0]) + int(hms.split(':')[1]) / 60.0 + float(hms.split(':')[2]) / 3600.0


# Translate RA,DEC sky coordinates to X,Y pixel values
def translate(ra, dec, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree):
    # Formulas based on IRAF implementation of rd2xy task
    det = cd11 * cd22 - cd12 * cd21
    if (det == 0.0):
        print >> stderr, '\n Error: Singular CD matrix'
        exit(-1)
    
    cdinv11 = cd22 / det
    cdinv12 = - cd12 / det
    cdinv21 = - cd21 / det
    cdinv22 = cd11 / det
    
    ra0 = radians( ra0 )
    dec0 = radians( dec0 )

    if degree == 'no':
        degra = hms2degree(ra)
        degdec = hms2degree(dec)
        degra = degra * 15.0
    else:
        degra = float(ra)
        degdec = float(dec)
        
    degra = radians(degra)
    degdec = radians(degdec)
    
    bottom = sin(degdec) * sin(dec0) + cos(degdec) * cos(dec0) * cos(degra - ra0)
    
    if (bottom == 0.0):
        print >> stderr, 'Error: Unreasonable RA/DEC range'
        exit(-1)

    xi = cos(degdec) * sin(degra - ra0) / bottom
    eta = (sin(degdec) * cos(dec0) - cos(degdec) * sin(dec0) * cos(degra - ra0)) / bottom
    xi = degrees(xi)
    eta = degrees(eta)
    
    x = cdinv11 * xi + cdinv12 * eta + crpix1
    y = cdinv21 * xi + cdinv22 * eta + crpix2
    
    print >> stdout, 'RA = %15s%s' %(ra, '\t'), '  DEC = %15s' %dec, '  X = %6.3f%s' %(x, '\t'), '  Y = %6.3f' %y
    
    return (ra, dec, x, y)
    
    
# sky2pix routine to proccess ra,dec pairs
def sky2pix(image, skyfile, degree = 'yes', outfile = None):
    # Read input image header keyword values
    ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22 = getHeader(image)  

    # Open input file
    try:
        ifile = open(skyfile, 'r')
    except:
        print >> stderr, 'Error: Not able to open input file ', skyfile, '. Exiting.'
        exit(-1)

    # Set the output file name
    if not outfile:
        if len(skyfile.rsplit('/')) > 1:
            outfile = join( skyfile.rsplit('/')[0], 'sky2pix.out' )
        else:
            outfile = 'sky2pix.out'

    # Open the output file
    try:    
        ofile = open(outfile, 'w')
    except:
        print >> stderr, 'Error: Not able to open output file ', outfile, '. Exiting.'
        exit(-1)
        
    ofile.write('# --------------------------------------------------------------------\n')
    ofile.write('#    RA            DEC                 X           Y      \n')
    ofile.write('# --------------------------------------------------------------------\n')
    
    # Process the input data and write to the output file
    while 1:
        line = ifile.readline()
        if not line:
            break
        if line[0] != '#':
            res = translate(line.split()[0], line.split()[1], ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree)
            ofile.write('%15s%15s%3s%6.3f%3s%6.3f%s' %(str(res[0]), str(res[1]), '\t', res[2], '\t', res[3], '\n'))

    # Close the input and the output files
    try:        
        ifile.close()
        ofile.close()            
    except:
        print >> stderr, 'Warning: Not able to close input and output file.'
        
    print >> stdout, '\n Results written to file - ', outfile


# Main function - doing some data validation before calling sky2pix method
# ========================================================================
def main(image, skyfile, degree = 'yes', outfile = None):
    if not exists(image.rsplit('[', 1)[0]):
        print >> stderr, 'Error: Image ', image, ' does not exist. Exiting.'
        exit(-1)

    if not exists(skyfile):
        print >> stderr, 'Error: Sky coordinate file ', skyfile, ' does not exist. Exiting.'
        exit(-1)  

    sky2pix(image, skyfile, degree, outfile)
    
    
# Entry point for SKY2PIX_SERIAL utility
# ======================================
if __name__ == '__main__':
    usage = "Usage: python %prog [options] image skyfile"
    description = "Description. Utility to convert RA/DEC sky coordinates to X/Y pixel image coordinates in serial mode."
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
    if len(args) != 2:
        parser.error("Incorrect number of arguments. Check help for more details.")

    print >> stdout, '\n Starting processing...'
    
    # Check verbosity
    if not options.verbose:
        output = StringIO()
        old_stdout = stdout
        stdout = output

    # Check if pyraf module is installed
    try:
        import pyfits
    except:
        print >> stderr, 'Error: Python module pyfits not found. Exiting.'
        exit(-1)
    
    main(args[0], args[1], options.degree, options.filename)
    
    # Reset verbosity
    if not options.verbose:
        stdout = old_stdout
    
    print >> stdout, '\n Process completed successfully.'