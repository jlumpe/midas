# testdb_210126 data

This directory is intended to contain data for the `testdb_210126` test database. Only the sqlite
file `testdb_210126.db` is included in version control and is used for tests of the SQLAlchemy
models and schema.

The other two files are used for the optional full query tests and need to be downloaded from the
MIDAS data repository (wherever that is):

* `testdb_210126.midas-signatures.gz`
* `testdb_210126_query_seqs.tar.gz` - extract to the `query_seqs/` subdirectory.
