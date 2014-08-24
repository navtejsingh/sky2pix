sky2pix
=======

Python routine to convert astronomical object sky coordinate (RA/DEC) to CCD pixel coordinates on multicore machines.

- Author:       Navtej Singh
- Contact:      reachnavtej@gmail.com
- Web Site:     http://astro.nuigalway.ie/staff/navtejs
- Organization: CfA@NUIG <http://astro.nuigalway.ie>

This routine was coded as part of research paper "Parallel astronomical data processing with Python: Recipes for multicore machines", published in Astronomy and Computing. Astro-ph link: http://arxiv.org/abs/1306.0573.

Thank you for downloading sky2pix code. It includes three different versions of the same code - sky2pix_serial.py for running code in serial mode,
sky2pix_multi.py for running on multicore/multiprocessor machines and sky2pix_pp.py uses parallel python module.

- Following requirements should be met to run the code.

    + A Python 2.4/2.5/2.6/2.7/3.0/3.1/3.2 distribution.

    + pyfits python module to handle FITS image files. Download from STScI's
      website (http://www.stsci.edu/resources/software_hardware/pyfits)

    + parallel python module. Can be downloaded from www.parallelpython.com      


- Execute the following commands to run the code

    + Serial Mode: 
                $python sky2pix_serial.py in.fits in_rd_deg.cat
                $python sky2pix_serial.py in.fits in_rd_hms.cat -d no 

    + Multicore Mode: 
                $python sky2pix_multi.py in.fits in_rd_deg.cat
                $python sky2pix_multi.py in.fits in_rd_hms.cat -d no    

    + Parallel python module: 
                $python sky2pix_pp.py in.fits in_rd_deg.cat
                $python sky2pix_pp.py in.fits in_rd_hms.cat -d no

    
- For available command line options:

    + $python sky2pix_serial.py --help
    + $python sky2pix_multi.py --help
    + $python sky2pix_pp.py --help