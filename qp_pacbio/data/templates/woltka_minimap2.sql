CREATE TEMP TABLE lengths AS
    SELECT *
        FROM read_csv('{{reference_len_map}}', header=false, delim='\t', columns = {"read_id": "VARCHAR", "length": "BIGINT"});

COPY (SELECT *
      FROM read_alignments('/dev/stdin', reference_lengths='lengths')
          where alignment_seq_identity(cigar, tag_nm, tag_md, 'blast') > 0.9
              and alignment_query_coverage(cigar) > 0.9
     ) TO '/dev/stdout' (FORMAT SAM, INCLUDE_HEADER false, COMPRESSION gzip);
