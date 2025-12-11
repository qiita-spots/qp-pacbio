COPY (SELECT *
      FROM read_alignments('/dev/stdin')
      WHERE alignment_query_coverage(cigar) >= 0.9
          AND alignment_seq_identity(cigar, tag_nm, tag_md, 'blast') >= 0.95
) TO '/dev/stdout' (FORMAT SAM, COMPRESSION gzip, INCLUDE_HEADER false);
