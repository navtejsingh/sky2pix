sky2pix
=======

Python routine to convert astronomical object sky coordinate (RA/DEC) to CCD pixel coordinates on multicore machines.

- Author:       Navtej Singh
- Contact:      reachnavtej@gmail.com
- Web Site:     http://astro.nuigalway.ie/staff/navtejs
- Organization: CfA@NUIG <http://astro.nuigalway.ie>

This routine was coded as part of research paper "Parallel astronomical data processing with Python: Recipes for multicore machines", published in Astronomy and Computing. Astro-ph link: http://arxiv.org/abs/1306.0573.

Thank you for downloading sky2pix code. It includes four different versions of the same code - sky2pix_serial.py for running code in serial mode,
sky2pix_multi.py for running on multicore/multiprocessor machines, sky2pix_pp.py uses parallel python module and sky2pix_mpi.py for running in parallel cluster mode.

- Following requirements should be met to run the code.

    + A Python 2.4/2.5/2.6/2.7/3.0/3.1/3.2 distribution.

    + pyfits python module to handle FITS image files. Download from STScI's
      website (http://www.stsci.edu/resources/software_hardware/pyfits)

    + parallel python module. Can be downloaded from www.parallelpython.com      

    + MPI library if running in parallel cluster mode. MPICH2 is recommended.
    
    + mpi4py module to communicate between python and MPI library

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
                
    + Parallel cluster Mode:
    
            $mpiexec -np ncpus python sky2pix.py in.fits in_rd_deg.cat
            $mpiexec -np ncpus python sky2pix.py in.fits in_rd_hms.cat -d no
			
		where ncpus is number of processors.
		
		It can also be executed using pbs_job if pbs/torque/maui is installed on the cluster.

    
- For available command line options:

    + $python sky2pix_serial.py --help
    + $python sky2pix_multi.py --help
    + $python sky2pix_pp.py --help
    + $mpiexec -np ncpus python sky2pix.py --help
