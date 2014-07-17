CREATE TABLE transient 
  (id SERIAL
  ,runcat INT NOT NULL
  ,band SMALLINT NOT NULL
  ,siglevel DOUBLE PRECISION DEFAULT 0
  ,v_int DOUBLE PRECISION NOT NULL
  ,eta_int DOUBLE PRECISION NOT NULL
  ,detection_level DOUBLE PRECISION DEFAULT 0
  ,trigger_xtrsrc INT NOT NULL
  ,transient_type SMALLINT NOT NULL
  -- NB previous_limits_image id can be constrained to not null when transient
  -- def goes dynamic and this becomes the 'newsource' table
  ,previous_limits_image INT DEFAULT NULL
  ,status INT DEFAULT 0
  ,t_start TIMESTAMP
{% ifdb postgresql %}
  ,PRIMARY KEY (id)
{% endifdb %}
  ,FOREIGN KEY (runcat) REFERENCES runningcatalog (id)
  ,FOREIGN KEY (band) REFERENCES frequencyband (id)
  ,FOREIGN KEY (trigger_xtrsrc) REFERENCES extractedsource (id)
  ,FOREIGN KEY (previous_limits_image) REFERENCES image (id)
);

{% ifdb postgresql %}
CREATE INDEX "transient_runcat" ON "transient" ("runcat");
CREATE INDEX "transient_band" ON "transient" ("band");
CREATE INDEX "transient_trigger_xtrsrc" ON "transient" ("trigger_xtrsrc");
{% endifdb %}
