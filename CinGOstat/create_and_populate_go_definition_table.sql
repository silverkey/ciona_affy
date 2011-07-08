CREATE TABLE go_definition(
  go_id VARCHAR(30) PRIMARY KEY,
  definition VARCHAR(256)
);

LOAD DATA LOCAL INFILE 'go_definition.txt' INTO TABLE go_definition;
