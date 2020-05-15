******************
Database Structure
******************


Concepts
========


Implementation
--------------

The database schema is defined using Object Relational Mapping (ORM)
through the `SQLAlchemy <https://www.sqlalchemy.org/>`__, meaning that
every database object corresponds to a Python class. The data itself is
stored in a standard relational database (SQLAlchemy supports many
DBMS's, in both versions of the application SQlite is used). SQLAlchemy
automatically translates commands in its Python API into SQL queries, so
there is no need to interact with the DBMS directly.


Versioning of database objects
------------------------------

It is anticipated that multiple versions of the database will be
released over time, so we need some way of tracking the changes of
specific database objects across versions. To that end most objects have
``key`` and ``key_version`` columns. ``key`` is a string that uniquely
identifies the object. The current convention is to use a
slash-separated hierarchical format, e.g.
``"refseq/assembly/GCF_000730045.1"`` for a genome. ``key_version`` is a
version number that is incremented when the object's data is updated in
a new release.


Database objects
================

.. py:currentmodule:: midas.database.base


Genome
------

ORM model base class: :class:`.Genome`

Corresponds to a single assembly (one or more contigs, but at least
partially assembled) from what should be a single sequencing run. The
same organism or strain may have several genome entries for it. Typically
this will correspond directly to a record in Genbank (assembly database).

This model simply stores the metadata of the genome, information on the
actual sequence data itself is stored in the ``sequence`` relationship.
A Genome entry has at most one sequence associated with it.

The data on this model should primarily pertain to the sample and
sequencing run itself. It would be updated if for example a better
assembly was produced from the original raw data, however more advanced
interpretation such as taxonomy assignments belong on an attached
`Genome annotations`_ object.


Sequence
--------

ORM model base class: :class:`.Sequence`

Stores metadata for genome's sequence data, if it is in the database.

Genomes may be present in the database with associated metadata as a
Genome object, but may not have a sequence stored. This contains
metadata for the sequence itself.


Genome set
----------

ORM model base class: :class:`.GenomeSet`

A collection of genomes along with additional annotations on each.

This will be used (among other things) to identify a set of genomes
a query will be run against. Like Genomes they can be distributed via
updates, in which case they should have a unique key and version number
to identify them.


Genome annotations
------------------

ORM model base class: :class:`.GenomeAnnotations`

Connects a genome with a genome set. Can
also carry additional annotations for the genome that are different
between sets. Mostly holds taxonomy information as that is frequently
a result of additional analysis on the sequence.


Kmer set collection
-------------------

ORM model base class: :class:`.KmerSetCollection`

A collection of `Kmer set`_\ s for a set of genomes calculated
with the same parameters.


Kmer set
--------

ORM model base class: :class:`.KmerSet`

Reference to a stored k-mer set for a genome.

The actual data can be retrieved in coordinate format with
:func:`.AbstractDatabase.load_kset_coords`.
