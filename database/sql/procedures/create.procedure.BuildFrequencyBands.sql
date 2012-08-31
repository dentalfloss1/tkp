--DROP PROCEDURE BuildFrequencyBands;

/**
 */
CREATE PROCEDURE BuildFrequencyBands()
BEGIN

  /* Some cataloged sources do not have spectral information
    (e.g. exoplanets) and are referred to freq 0
   */
  INSERT INTO frequencyband
    (id
    ,freq_central
    ,freq_low
    ,freq_high
    ) 
  VALUES 
    (0
    ,NULL
    ,NULL
    ,NULL
    )
  ;

  INSERT INTO frequencyband
    (freq_central
    ,freq_low
    ,freq_high
    ) 
  VALUES 
    (31000000
    ,30000000 
    ,32000000
    )
    ,
    (37000000
    ,36000000 
    ,38000000
    )
    ,
    (43000000
    ,42000000 
    ,44000000
    )
    ,
    (49000000
    ,48000000 
    ,50000000
    )
    ,
    (54000000
    ,53000000 
    ,55000000
    )
    ,
    (60000000
    ,59000000 
    ,61000000
    )
    ,
    (66000000
    ,65000000 
    ,67000000
    )
    ,
    (74000000
    ,73000000 
    ,75000000
    )
    ,
    (120000000
    ,120000000 - 350000 / 2
    ,120000000 + 350000 / 2
    )
    ,
    (130000000
    ,130000000 - 450000 / 2
    ,130000000 + 450000 / 2
    )
    ,
    (140000000
    ,140000000 - 550000 / 2
    ,140000000 + 550000 / 2
    )
    ,
    (150000000
    ,150000000 - 700000 / 2
    ,150000000 + 700000 / 2
    )
    ,
    (160000000
    ,160000000 - 850000 / 2
    ,160000000 + 850000 / 2
    )
    ,
    (170000000
    ,170000000 - 1100000 / 2
    ,170000000 + 1100000 / 2
    )
    ,
    (325000000
    ,325000000 - 10000000 / 2
    ,325000000 + 10000000 / 2
    )
    ,
    (352000000
    ,352000000 - 20000000 / 2
    ,352000000 + 20000000 / 2
    )
    ,
    (640000000
    ,640000000 - 100000000 / 2
    ,640000000 + 100000000 / 2
    )
    ,
    (850000000
    ,850000000 - 100000000 / 2
    ,850000000 + 100000000 / 2
    )
    ,
    (1400000000
    ,1400000000 - 260000000 / 2
    ,1400000000 + 260000000 / 2
    )
    ,
    (2300000000
    ,2300000000 - 250000000 / 2
    ,2300000000 + 250000000 / 2
    )
    ,
    (4800000000
    ,4800000000 - 250000000 / 2
    ,4800000000 + 250000000 / 2
    )
    ,
    (8500000000
    ,8500000000 - 250000000 / 2
    ,8500000000 + 250000000 / 2
    )
  ;

END;
