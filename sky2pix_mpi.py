"""
   Routine to translate sky coordinates to x,y pixel values in multiprocessing mpi mode.
   Called from sky2pix.py.
"""

# Load python modules to be used in the routine
import os, sys
from StringIO import StringIO
from math import cos, sin, radians, degrees


# Check if multiprocesssing module is present
try:
    from mpi4py import MPI
except ImportError:
    raise ImportError("Module mpi4py not found.")


# Create chunks of data to be distributed to multiple processors based on file size
def getchunks(infile, n_cpus, scheduler = 'guided'):
    # Divide input data based on scheduler type
    if scheduler == 'static':
        size = os.path.getsize(infile) / n_cpus
    else:
        size = os.path.getsize(infile) / (n_cpus * 20)

    # Open input file    
    try:    
        ifile = open( infile )
    except:
        print >> sys.stderr, 'Error: Not able to open ', infile, '. Exiting.'
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
        print >> sys.stderr, 'Warning: Error closing the file ', ifile
    

# Get header keywords
def getHeader(image):
    print >> sys.stdout, '\n Getting image header keywords...'

    # Check if pyraf module is installed
    try:
        import pyfits
    except:
        print >> sys.stderr, 'Error: Python module pyfits not found. Exiting.'
        exit(-1)
    
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

    # Return image header keyword values
    return ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22


# Convert RA,DEC from hour:min:sec/degree:arcmin:arcsec to degree
def hms2degree(hms):
    return int(hms.split(':')[0]) + int(hms.split(':')[1]) / 60.0 + float(hms.split(':')[2]) / 3600.0


# Translate RA,DEC sky coordinates to X,Y pixel values
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
        print >> sys.stderr, 'Error: Unreasonable RA/DEC range'
        sys.exit(-1)

    xi = cos(degdec) * sin(degra - ra0) / bottom
    eta = (sin(degdec) * cos(dec0) - cos(degdec) * sin(dec0) * cos(degra - ra0)) / bottom
    xi = degrees(xi)
    eta = degrees(eta)
    
    x = cdinv11 * xi + cdinv12 * eta + crpix1
    y = cdinv21 * xi + cdinv22 * eta + crpix2
    
    print >> sys.stdout, 'RA = %15s%s' %(ra, '\t'), '  DEC = %15s' %dec, '  X = %6.3f%s' %(x, '\t'), '  Y = %6.3f' %y
    
    return (ra, dec, x, y)


# MPI worker function
def worker(indata):
    # Unpack input python list
    chunk0, chunk1, infile, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree = indata
    
    # Open input file
    try:
        ifile = open(infile, 'r')
    except:
        print >> sys.stderr, 'Error: Not able to open input file ', infile, '. Exiting.'
        sys.exit(-1)
        
    ifile.seek(chunk0)
    
    # Process data and save in list
    result = []
    for line in ifile.read(chunk1).splitlines():
    	if line[0] != '#':
    	    result.append(translate(line.split()[0], line.split()[1], ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree))

    # Return python result list        
    return result


# sky2pix to process input ra,dec pairs
def sky2pix(image, infile, degree = 'yes', outfile = None, verbose = True):

    # Start mpi processing
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    n_cpus = comm.Get_size()

    # Divide the data in chunks equal to number of processors
    if rank == 0:
        print >> sys.stdout, '\n Starting processing...'
        
        if not os.path.exists(image.rsplit('[', 1)[0]):
            print >> sys.stderr, 'Error: Image ', image, ' does not exist. Exiting.'
            sys.exit(-1)

        if not os.path.exists(infile):
            print >> sys.stderr, 'Error: Pixel file ', infile, ' does not exist. Exiting.'
            sys.exit(-1)

        # Read image headers
        ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22 = getHeader(image)	


        # Divide data in smaller chunks for better performance
        chunks = getchunks(infile, n_cpus, 'static')
    
        data = []
        for chunk in chunks:
            data.append((chunk[0], chunk[1], infile, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree))
    else:
        data = None

    # Check verbosity
    if not verbose:
        output = StringIO()
        old_stdout = sys.stdout
        sys.stdout = output

    # Scatter the data to available processors and gather the results
    indata = comm.scatter(data, root = 0)
    outdata = worker(indata)
    results = comm.gather(outdata, root = 0)

    # Write results to the output file
    if rank == 0:
        # Write headers to the output file
        if outfile:
            outfile = outfile
        else:
            outfile = os.path.join(infile.rsplit('/')[0], 'sky2pix.out')

        # Open output file for writing    
        try:
            ofile = open(outfile, 'w')
        except:
            print >> sys.stderr, 'Error: Not able to open file ', outfile, ' to write. Exiting.'
            sys.exit(-1)

        ofile.write('# --------------------------------------------------------------\n')
        ofile.write('#    X		Y		    RA	      	        DEC  \n')
        ofile.write('# --------------------------------------------------------------\n')
        
        for result in results:
            for value in result:
                ofile.write('%15s%15s%3s%6.3f%3s%6.3f%s' %(str(value[0]), str(value[1]), '\t', value[2], '\t', value[3], '\n'))

        # Close output file        
        try:        
            ofile.close()
        except:
            print >> sys.stderr, 'Warning: Not able to open file ', outfile
            
        
        print >> sys.stdout, '\n Results written to - ', outfile


    # Reset verbosity
    if not verbose:
        sys.stdout = old_stdout
    
    if rank == 0:
        print >> sys.stdout, '\n Process completed successfully.'
