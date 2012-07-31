#                                             LOFAR TRANSIENT DETECTION PIPELINE
#
#                                                         Source extraction node
#                                                                Evert Rol, 2011
#                                                          evert.astro@gmail.com
# ------------------------------------------------------------------------------

from __future__ import with_statement

import os
import sys
from threading import Thread
from contextlib import closing

from lofar.parameterset import parameterset
from lofarpipe.support.lofarnode import LOFARnodeTCP
from lofarpipe.support.utilities import log_time
from lofarpipe.support.lofarexceptions import PipelineException

def store_to_mongodb(filename, hostname, port, db, logger):
    logger.info(
        "Storing %s to MongoDB database %s on %s port %d" %
        (filename, db, hostname, port)
    )
    try:
        import pymongo
        import gridfs
    except ImportError:
        logger.warn("Could not import MongoDB modules")
        return

    try:
        connection = pymongo.Connection(host=hostname, port=port)
        gfs = gridfs.GridFS(connection[db])
        new_file = gfs.new_file(filename=filename)
        with open(filename, "r") as f:
            new_file.write(f)
        new_file.close()
        connection.close()
    except Exception, e:
        logger.warn("Could not store image to MongoDB: %s" % (str(e)))

class source_extraction(LOFARnodeTCP):
    """
    Extract sources from a FITS image
    """

    def run(
        self, image, dataset_id, parset,
        store_images, mongo_host, mongo_port, mongo_db,
        tkpconfigdir=None
    ):
        """
        Args:

            - image: FITS filename

            - dataset_id: dataset to which image belongs

            - parset: parameter set *filename* containg at least the
                  detection threshold and the source association
                  radius, the last one in units of the de Ruiter
                  radius.

            - storage_images: bool. Store images to MongoDB database if True.

            - mongo_host/port/db: details of MongoDB to use if store_images is
              True.
        """
        if tkpconfigdir:   # allow nodes to pick up the TKPCONFIGDIR
            os.environ['TKPCONFIGDIR'] = tkpconfigdir
        import tkp
        from tkp.config import config
        from tkp.database import DataBase, DataSet
        from tkp.utility.accessors import FITSImage
        from tkp.utility.accessors import dbimage_from_accessor
        from tkp.utility.accessors import sourcefinder_image_from_accessor

        if not os.path.exists(image):
            raise PipelineException("Image '%s' not found" % image)

        with log_time(self.logger):
            with closing(DataBase()) as database:
                seconfig = config['source_extraction']
                parset = parameterset(parset)
                dataset = DataSet(id=dataset_id, database=database)
                fitsimage = FITSImage(image)
                db_image = dbimage_from_accessor(dataset=dataset,
                                                 image=fitsimage)
                self.logger.info("Detecting sources in %s at %f level",
                                 image, parset.getFloat('detection.threshold'))
                
                data_image = sourcefinder_image_from_accessor(fitsimage)
                if store_images:
                    upload_thread = Thread(
                        target=store_to_mongodb,
                        args=(image, mongo_host, mongo_port, mongo_db, self.logger)
                    )
                    upload_thread.start()

                seconfig['back_sizex'] = parset.getInt('backsize.x',
                                                       seconfig['back_sizex'])
                seconfig['back_sizey'] = parset.getInt('backsize.y',
                                                       seconfig['back_sizey'])
                seconfig['margin'] = parset.getFloat('margin',
                                                   seconfig['margin'])
                seconfig['deblend'] = parset.getBool('deblend',
                                                   seconfig['deblend'])
                seconfig['deblend_nthresh'] = parset.getInt('deblend_nthresh',
                                                   seconfig['deblend_nthresh'])
                
                
                self.logger.info("Employing margin: %f, deblend: %s, deblend_nthresh:%d",
                                 seconfig['margin'], 
                                 seconfig['deblend'],
                                 seconfig['deblend_nthresh']
                                 )  
                
                det = parset.getFloat('detection.threshold',
                                      seconfig['detection_threshold'])
                anl = parset.getFloat('analysis.threshold',
                                      seconfig['analysis_threshold'])
                
                ##Finally, do some work!
                results = data_image.extract(det=det, anl=anl)
                
                self.logger.info("Detected %d sources", len(results))
                self.logger.info("Saving extracted sources to database")
                tuple_results = [result.serialize() for result in results]
                db_image.insert_extracted_sources(tuple_results)
                deRuiter_r = (parset.getFloat('association.radius') *
                              config['source_association']['deruiter_radius'])
                self.logger.info("Associate newly extracted sources with existing ones")
                db_image.associate_extracted_sources(deRuiter_r=deRuiter_r)
                #self.logger.info("Update monitoring list for already found sources")
                #db_image.match_monitoringlist(assoc_r=deRuiter_r, mindistance=30)
                self.outputs['image_id'] = db_image.id

                if store_images:
                    upload_thread.join()
        return 0

if __name__ == "__main__":
    #   If invoked directly, parse command line arguments for logger information
    #                        and pass the rest to the run() method defined above
    # --------------------------------------------------------------------------
    jobid, jobhost, jobport = sys.argv[1:4]
    sys.exit(source_extraction(jobid, jobhost, jobport).run_with_stored_arguments())
