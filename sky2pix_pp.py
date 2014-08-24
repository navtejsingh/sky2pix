#!/usr/bin/env python

"""
    ------------------------------------------------------------------------
    Routine to translate RA and DEC sky coordinates to X, Y pixel values 
    using parallel python.

    Usage: python sky2pix_pp.py [options] image skyfile

    Input:
        image: Input image with basic header keywords
        skyfile: Input file with (ra, dec) pixel value in column 1 and 2
      
      
    [Options]:
        --help: help
        --version: program version
        --verbose: show result messages
        --quiet: don't show result messages
        --filename: output file name (default is pix2sky.dat)
        --hour: (ra, dec) in hms? (default is yes)
        --ncpus: Number of processors to use [default is maximum number of cores]
        --scheduler: Type of scheduler (static, guided, dynamic) [default is guided]
        
        
    Output:
        sky2pix.out: Output file with (ra, dec) and corresponding (x, y)
        pixel values
        
    Author:
        Navtej Singh

    Organization:
        Centre for Astronomy, National University of Ireland, Galway, Ireland

    Version:
        24 October 2012     1.0     Original version
    ------------------------------------------------------------------------            
"""


# Load python modules to be used in the routine
import sys, math
from os.path import getsize, join, exists
from StringIO import StringIO
from optparse import OptionParser


# Check if parallel python module is available
try:
    import pp
except:
    print >> sys.stderr, 'Error: Parallel python module not found. Use serial version of this routine. Exiting.'
    sys.exit(-1)


# Create chunks of data to be distributed to multiple processors 
# based on the file size. Arbitrary factor of 20 was chosen
def getchunks(infile, n_cpus, scheduler = 'guided'):
    # Divide the input data based on the scheduler type
    if scheduler == 'static':
        size = getsize(infile) / n_cpus
    else:
        size = getsize(infile) / ( n_cpus * 20 )

    # Open the input file    
    try:    
        ifile = open(infile)
    except:
        print >> sys.stderr, 'Error: Not able to open ', infile, '. Exiting.'
        sys.exit(-1)

    # Create chunk of data to be distributed to the nodes    
    while 1:
        start = ifile.tell()
        ifile.seek(size, 1)
        s = ifile.readline()
        yield start, ifile.tell() - start
        if not s:
            break

    # Close the input file    
    try:    
        ifile.close()
    except:
        print >> sys.stderr, 'Warning: Error closing the file ', ifile
    

# Get the image header keywords
# =============================
def getHeader(image):
    print >> sys.stdout, '\n Getting image header keywords...'
    
    # Input can be a single FITS image, multi-extension image or
    # multi-extension image with particular extension
    if len(image.split('[', 1)) > 1:
        ext = image.split('[', 1)[1].replace(']', '')
        image = image.split('[', 1)[0]
    else:
        ext = ''
    
    # Open Header Unit List (HDU) to read header keywords
    try:
        hdulist = pyfits.open(image)
    except:
        print >> sys.stderr, 'Error: Not able to read FITS header. Exiting.'
        sys.exit(-1)

    hdulist.close()

    # Get header parameters - checking number of extensions and using
    # the 1st extension in case of multi extension FITS image
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

    # Return image header keyword values to the calling method
    return ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22


# Convert RA from hour:min:sec to degree
# ======================================
def hms2degree(hms):
    return int(hms.split(':')[0]) + int(hms.split(':')[1]) / 60.0 + float(hms.split(':')[2]) / 3600.0


# Translate RA,DEC sky coordinates to X,Y pixel values
# ====================================================
def translate(ra, dec, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree):
    # Formulas based on IRAF implementation of rd2xy task
    det = cd11 * cd22 - cd12 * cd21
    if (det == 0.0):
        print >> sys.stderr, '\n Error: Singular CD matrix'
        sys.exit(-1)
    
    cdinv11 = cd22 / det
    cdinv12 = - cd12 / det
    cdinv21 = - cd21 / det
    cdinv22 = cd11 / det
    
    ra0 = math.radians(ra0)
    dec0 = math.radians(dec0)

    if degree == 'no':
        degra = hms2degree(ra)
        degdec = hms2degree(dec)
        degra = degra * 15.0
    else:
        degra = float(ra)
        degdec = float(dec)
        
    degra = math.radians(degra)
    degdec = math.radians(degdec)
    
    bottom = math.sin(degdec) * math.sin(dec0) + math.cos(degdec) * math.cos(dec0) * math.cos(degra - ra0)
    
    if (bottom == 0.0):
        print >> sys.stderr, 'Error: Unreasonable RA/DEC range'
        sys.exit(-1)

    xi = math.cos(degdec) * math.sin(degra - ra0) / bottom
    eta = (math.sin(degdec) * math.cos(dec0) - math.cos(degdec) * math.sin(dec0) * math.cos(degra - ra0)) / bottom
    xi = math.degrees(xi)
    eta = math.degrees(eta)
    
    x = cdinv11 * xi + cdinv12 * eta + crpix1
    y = cdinv21 * xi + cdinv22 * eta + crpix2
    
    print >> sys.stdout, 'RA = %15s%s' %(ra, '\t'), '  DEC = %15s' %dec, '  X = %6.3f%s' %(x, '\t'), '  Y = %6.3f' %y
    
    return (ra, dec, x, y)


# Parallel python worker function
# ===============================
def worker(indata):
    # Unpack input python list
    chunk0, chunk1, infile, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree = indata

    # Open the input file
    try:
        ifile = open(infile, 'r')
    except:
        print >> sys.stderr, 'Error: Not able to open the input file ', infile, '. Exiting.'
        sys.exit(-1)

    # Process input records    
    ifile.seek(chunk0)
    result = []
    for line in ifile.read(chunk1).splitlines():
        if line[0] != '#':
            result.append(translate(line.split()[0], line.split()[1], ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree))

    # Return results to calling method        
    return result

       
# sky2pix routine to process ra, dec pairs
# ========================================
def sky2pix(image, infile, degree = 'yes', outfile = None, ncpus = None, scheduler = 'guided'):
    # Read image header keyword values
    ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22 = getHeader(image)  

    # Allocate number of worker proceses equal to usre input or default to 
    # the maximum number of processors/cores on the machine
    ppservers = ()
    if ncpus:
        # Creates jobserver with ncpus workers
        job_server = pp.Server(int(ncpus), ppservers=ppservers)
    else:
        # Creates jobserver with automatically detected number of workers
        job_server = pp.Server(ppservers = ppservers)
    
    # Divide data in smaller chunks for better performance based on scheduler type
    chunks = getchunks(infile, job_server.get_ncpus(), scheduler)

    # Start the worker processes in parallel
    jobs = []
    for value in chunks:
            indata = (value[0], value[1], infile, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree)
            jobs.append(job_server.submit(worker, (indata,), (translate,hms2degree,), ("math","sys",)))

    # Append the output results
    results = []
    for job in jobs:
        results.append(job())
    

    # Set the output file name
    if not outfile:
        if len(infile.rsplit('/')) > 1:
            outfile = join(infile.rsplit('/')[0], 'sky2pix.out')
        else:
            outfile = 'sky2pix.out'  
            
    # Open the output file for writing    
    try:    
        ofile = open(outfile, 'w')
    except:
        print >> sys.stderr, 'Error: Not able to open output file ', outfile, '. Exiting.'
        sys.exit(-1)

    # Write data header    
    ofile.write('# --------------------------------------------------------------\n')
    ofile.write('#    X     Y           RA                  DEC  \n')
    ofile.write('# --------------------------------------------------------------\n')

    # Write data to the output file
    for result in results:
        for value in result:
            ofile.write('%15s%15s%3s%6.3f%3s%6.3f%s' %(str(value[0]), str(value[1]), '\t', value[2], '\t', value[3], '\n'))   
    
    # Close the output file     
    try:        
        ofile.close()
    except:
        print >> sys.stderr, 'Warning: Not able to close output file ', outfile
        
    print >> sys.stdout, '\n Results written to - ', outfile



# Main function - doing some data validation before calling sky2pix method
# ========================================================================
def main(image, skyfile, degree = 'yes', outfile = None, ncpus = None, scheduler = 'guided'):
    # Check if the image exists
    if not exists(image.rsplit('[', 1)[0]):
        print >> sys.stderr, 'Error: Image ', image, ' does not exist. Exiting.'
        sys.exit(-1)

    # Check if the input sky value file exists 
    if not exists(skyfile):
        print >> sys.stderr, 'Error: Pixel file ', skyfile, ' does not exist. Exiting.'
        sys.exit(-1)  

    # Execute the main method
    sky2pix(image, skyfile, degree, outfile, ncpus, scheduler)
    
    
# Entry point for SKY2PIX_MULTI utility
# =====================================
if __name__ == '__main__':
    usage = "Usage: python %prog [options] image skyfile"
    description = "Description. Utility to convert RA/DEC sky coordinates to X/Y pixel image coordinates in multiprocessing mode."
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
    parser.add_option("-n", "--ncpus", dest = "ncpus", metavar="NCPUS",
                    action="store", help = "number of cpus (cores) for processing"
                    )
    parser.add_option("-s", "--scheduler", dest = "scheduler", metavar="SCHEDULER",
                    action="store", help = "scheduler for multiprocessing [default is guided]",
                    choices=['guided', 'static'], default = 'guided'
                    )
    (options, args) = parser.parse_args()
    
    # Check for number of input arguments
    if len(args) != 2:
        parser.error("Incorrect number of arguments. Check help for more details.")

    print >> sys.stdout, '\n Starting processing...'
    
    # Check verbosity
    if not options.verbose:
        output = StringIO()
        old_stdout = sys.stdout
        sys.stdout = output

    # Check if python pyfits module is available
    try:
        import pyfits
    except:
        print >> sys.stderr, 'Error: Python module pyfits not found. Exiting.'
        sys.exit(-1)
    
    main(args[0], args[1], options.degree, options.filename, options.ncpus, options.scheduler)
    
    # Reset verbosity
    if not options.verbose:
        sys.stdout = old_stdout
    
    print >> sys.stdout, '\n Process completed successfully.'