Collapse the records of a SAM file

Reads that have exactly the same coordinates are collapsed into a single entry. A new tag that reports the number of records that have been collapsed is added to each new collapsed record. If exists, the quality field is set to "*" and the name of the first record is used for the new collapsed record.
