CREATE TABLE dataset
  (id INT AUTO_INCREMENT
  ,rerun INT NOT NULL DEFAULT '0'
  ,"type" TINYINT NOT NULL
  ,process_ts TIMESTAMP NOT NULL
  ,detection_threshold DOUBLE NOT NULL
  ,analysis_threshold DOUBLE NOT NULL
  ,assoc_radius DOUBLE NOT NULL
  ,backsize_x SMALLINT NOT NULL
  ,backsize_y SMALLINT NOT NULL
  ,margin_width DOUBLE NOT NULL
  ,inname VARCHAR(64) NOT NULL
  ,description VARCHAR(100) DEFAULT NULL
  ,node TINYINT NOT NULL DEFAULT %NODE%
  ,nodes TINYINT NOT NULL DEFAULT %NODES%
  ,PRIMARY KEY (id)
  )
;

