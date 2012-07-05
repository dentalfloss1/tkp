--DROP FUNCTION insertImage;

/**
 * This function inserts a row in the image table,
 * and returns the value of the id under which it is known.
 * If the dataset name (dsinname) already exists, a new row is added
 * and the rerun value is incremented by 1. If not, it is set
 * to its default value 0.
 *
 * Note I: To be able to create a function that modifies data 
 * (by insertion) we have to set the global bin log var:
 * mysql> SET GLOBAL log_bin_trust_function_creators = 1;
 *
 * Note II: The params in comment should be specified soon.
 * This means this function inserts deafult values so long.
 *
 */
CREATE FUNCTION insertImage(ids_id INT
                           ,itau INT
                           ,itau_time DOUBLE
                           ,ifreq_eff DOUBLE
                           ,ifreq_bw DOUBLE
                           ,itaustart_ts TIMESTAMP
                           ,ibeam_maj DOUBLE
                           ,ibeam_min DOUBLE
                           ,ibeam_pa DOUBLE
                           ,iurl VARCHAR(1024)
                           ) RETURNS INT
BEGIN

  DECLARE iimageid INT;
  DECLARE oimageid INT;
  /*DECLARE iseq_nr INT;*/
  DECLARE iband INT;

  SET iband = getBand(ifreq_eff);
  
  /*
  SELECT MAX(seq_nr)
    INTO iseq_nr
    FROM images
   WHERE ds_id = ids_id
     AND tau = itau
     AND band = iband
  ;

  IF iseq_nr IS NULL THEN
    SET iseq_nr = 1;
  ELSE
    SET iseq_nr = iseq_nr + 1;
  END IF;
  */

  SELECT NEXT VALUE FOR seq_images INTO iimageid;

  INSERT INTO images
    (imageid
    ,ds_id
    ,tau
    ,band
    ,tau_time
    ,freq_eff
    ,freq_bw
    ,taustart_ts
    ,bsmaj
    ,bsmin
    ,bpa
    ,url
    ) 
  VALUES
    (iimageid
    ,ids_id
    ,itau
    ,iband
    ,itau_time
    ,ifreq_eff
    ,ifreq_bw
    ,itaustart_ts
    ,ibeam_maj 
    ,ibeam_min 
    ,ibeam_pa 
    ,iurl
    )
  ;

  SET oimageid = iimageid;
  RETURN oimageid;

END;
