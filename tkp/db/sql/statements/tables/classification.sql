/**
 * TODO: evert's tabel, how should we handle this? what is the ID?
 */

CREATE TABLE classification
  (id SERIAL
	,transient_id INTEGER NOT NULL
	,classification VARCHAR(256)
	,weight DOUBLE PRECISION NOT NULL DEFAULT 0

{% ifdb postgresql %}
  ,PRIMARY KEY (id)
{% endifdb %}
);
