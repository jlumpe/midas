File Formats
============

There are several file formats used to store and distribute MIDAS data
outside of an actual database. The purpose of these different formats is
to store the information in a somewhat modular and platform-independent
way. For example, genome data is stored separately from thresholds, and genome
data is just a bunch of JSON files in a zip archive. This makes the data
relatively transparent and does not rely on any particular database system.


Archive files
-------------

These have the ``.midas-archive`` extension and store one or more genome
sets along with their genomes and annotatios. They are actually zip files
containing JSON files with a specific directory structure. Use the
:class:`midas.archive.DatabaseArchive` class to read these and import them into
an existing database.


Thresholds
----------

.. py:currentmodule:: midas.database.base

These are just serialized Python dictionaries in
`pickle <https://docs.python.org/3/library/pickle.html>`_ format. They have
a ``'species'`` and/or ``'genus'`` item, which is another dictionary that maps
species/genus keys to threshold values. The keys are taken from the
``tax_genus`` and ``tax_species`` fields of the :class:`.GenomeAnnotations`
model (both strings). Genus dictionaries are keyed by genus string only,
species dictionaries are keyed by ``(genus, species)`` tuples.

Threshold files contain no metadata like key/version fields and are not imported
into a database like the contents of archive files (they are used as-is).
Therefore, they are only distinguished by file name.