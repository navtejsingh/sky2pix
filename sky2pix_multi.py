#!/usr/bin/env python

"""
    ------------------------------------------------------------------------
    Routine to translate RA and DEC sky coordinates to X, Y pixel values in
    multiprocessing mode

    Usage: python sky2pix_multi.py [options] image skyfile

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
        --mode: multiprocessing mode (pool, process) [default is process]
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
        10 January 2012     1.0     Original version
    ------------------------------------------------------------------------            
"""


# Load python modules to be used in the routine
from os.path import getsize, join, exists
from sys import stderr, stdout, exit
from StringIO import StringIO
from math import cos, sin, radians, degrees
from optparse import OptionParser


# Check if multiprocesssing module is present
try:
    import multiprocessing as mp
    from multiprocessing import Process, Queue
except:
    print >> stderr, 'Error: Python multiprocessing module not found. Use serial version of this routine. Exiting.'
    exit(-1)


# Create chunks of data to be distributed to multiple processors 
# based  on the file size. Arbitrary factor of 20 was chosen.
def getchunks(infile, n_cpus, scheduler = 'guided'):
    # Divide input data based on scheduler type
    if scheduler == 'static':
        size = getsize(infile) / n_cpus
    else:
        size = getsize(infile) / (n_cpus * 20)

    # Open input file    
    try:    
        ifile = open(infile)
    except:
        print >> stderr, 'Error: Not able to open ', infile, '. Exiting.'
        exit(-1)

    # Create chunk of data to be distributed to nodes    
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
        print >> stderr, 'Warning: Error closing the file ', ifile
    

# Get header keywords
# ===================
def getHeader(image):
    print >> stdout, '\n Getting image header keywords...'
    
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
        print >> stderr, 'Error: Not able to read FITS header. Exiting.'
        exit(-1)

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
        print >> stderr, '\n Error: Singular CD matrix'
        exit(-1)
    
    cdinv11 = cd22 / det
    cdinv12 = - cd12 / det
    cdinv21 = - cd21 / det
    cdinv22 = cd11 / det
    
    ra0 = radians(ra0)
    dec0 = radians(dec0)

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


# Multicore process worker function  
# ================================= 
def process_worker(s_q, r_q, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree):
    for value in iter(s_q.get, 'STOP'):
        # Open input data file
        try:
            ifile = open(value[0], 'r')
        except:
            print >> stderr, 'Error: Not able to open input file ', value[0], '. Exiting.'
            exit(-1)

        # Process input record and write to output python list    
        ifile.seek(value[1])
        result = []
        for line in ifile.read(value[2]).splitlines():
            if line[0] != '#':
                result.append(translate(line.split()[0], line.split()[1], ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree))

        # Put results in receive queue        
        r_q.put(result)


# Multicore pool worker function
# ==============================
def pool_worker(indata):
    # Unpack input python list
    chunk0, chunk1, infile, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree = indata

    try:
        ifile = open(infile, 'r')
    except:
        print >> stderr, 'Error: Not able to open input file ', infile, '. Exiting.'
        exit(-1)

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
def sky2pix(image, infile, degree = 'yes', outfile = None, mode = 'process', ncpus = None, scheduler = 'guided'):
    # Read image header keyword values
    ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22 = getHeader(image)  

    # Get number of processors (cores) on the machine. In case of Intel processors with
    # hyperthreading, total number of processors will be equal to number of cores * number of threads/core
    if ncpus:
        n_cpus = int(ncpus)
    else:
        n_cpus = mp.cpu_count()
    
    
    # Divide data in smaller chunks for better performance
    chunks = getchunks(infile, n_cpus, scheduler)
    
    
    # Based on multiprocessing mode, follow seperate path
    if mode == 'process':
        send_q = Queue()
        recv_q = Queue()
    
        # Start number of processes equal to number of cores
        for i in range(n_cpus):
            Process(target = process_worker, args = (send_q, recv_q, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree)).start()
        
        # Put data in the send queue
        cnt = 0
        for chunk in chunks:
            cnt += 1
            send_q.put((infile, chunk[0], chunk[1]))
        
        # Set the output file name
        if not outfile:
            if len(infile.rsplit('/')) > 1:
                outfile = join(infile.rsplit('/')[0], 'sky2pix.out')
            else:
                outfile = 'sky2pix.out'   

        # Open the output file
        try:
            ofile = open(outfile, 'w')
        except:
            print >> stderr, 'Error: Not able to open the outfile file ', outfile, '. Exiting.'
            exit(-1)

        # Write headers to the output file     
        ofile.write('# ----------------------------------------------------------------\n')
        ofile.write('#    RA            DEC              X              Y  \n')
        ofile.write('# ----------------------------------------------------------------\n')

        # Write results to output file
        for i in range(cnt):
            res = recv_q.get()
            for value in res:
                ofile.write('%15s%15s%3s%6.3f%3s%6.3f%s' %(str(value[0]), str(value[1]), '\t', value[2], '\t', value[3], '\n'))

        # Close the output file        
        try:        
            ofile.close()
        except:
            print >> stderr, 'Warning: Not able to close output file ', outfile
        
        # Stop all the running processes
        for i in range(n_cpus):
            send_q.put('STOP')
            
        print >> stdout, '\n Results written to - ', outfile
    else:
        # Create python list to pass pool processes
        indata = []
        for chunk in chunks:
            indata.append((chunk[0], chunk[1], infile, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree))

        # Set the output file name
        if not outfile:
            if len(infile.rsplit('/')) > 1:
                outfile = join(infile.rsplit('/')[0], 'sky2pix.out')
            else:
                outfile = 'sky2pix.out'  

        # Open output file for writing    
        try:    
            ofile = open(outfile, 'w')
        except:
            print >> stderr, 'Error: Not able to open output file ', outfile, '. Exiting.'
            exit(-1)

        # Write data header    
        ofile.write('# --------------------------------------------------------------\n')
        ofile.write('#    X     Y           RA                  DEC  \n')
        ofile.write('# --------------------------------------------------------------\n')

        
        # To prevent memory overflow, process and write data in chunks
        for i in xrange(0, len(indata), n_cpus):
            tmpdata = indata[i : i + n_cpus]
            pool = mp.Pool(n_cpus)
            result = pool.map(pool_worker, tmpdata)
    
            for res in result:
                for value in res:
                    ofile.write('%15s%15s%3s%6.3f%3s%6.3f%s' %(str(value[0]), str(value[1]), '\t', value[2], '\t', value[3], '\n'))

        # Close the output file
        try:        
            ofile.close()
        except:
            print >> stderr, 'Warning: Not able to close output file ', outfile
        
        print >> stdout, '\n Results written to - ', outfile


# Main function - doing some data validation before calling sky2pix method
def main(image, skyfile, degree = 'yes', outfile = None, mode = 'process', ncpus = None, scheduler = 'guided'):
    # Check of the image exists
    if not exists(image.rsplit('[', 1)[0]):
        print >> stderr, 'Error: Image ', image, ' does not exist. Exiting.'
        exit(-1)

    # Check if the input sky value file exists
    if not exists(skyfile):
        print >> stderr, 'Error: Sky value file ', skyfile, ' does not exist. Exiting.'
        exit(-1)  

    # Execute the method
    sky2pix(image, skyfile, degree, outfile, mode, ncpus, scheduler)
    
    
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
    parser.add_option("-m", "--mode", dest = "mode", metavar="MODE",
                    action="store", help = "multiprocessing mode (pool or process) [default is process]",
                    choices=['process', 'pool'], default = 'process'
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
    
    main(args[0], args[1], options.degree, options.filename, options.mode, options.ncpus, options.scheduler)
    
    # Reset verbosity
    if not options.verbose:
        stdout = old_stdout
    
    print >> stdout, '\n Process completed successfully.'