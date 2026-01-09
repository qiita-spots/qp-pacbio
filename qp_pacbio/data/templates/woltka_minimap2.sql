CREATE TEMP TABLE lengths AS
    SELECT *
        FROM read_csv('{{reference_len_map}}', header=false, delim='\t', columns = {"read_id": "VARCHAR", "length": "BIGINT"});

COPY (SELECT *
      FROM read_alignments('/dev/stdin', reference_lengths='lengths')
          where alignment_seq_identity(cigar, tag_nm, tag_md, 'blast') > {{identity}}
              and alignment_query_coverage(cigar) > {{coverage}}
     ) TO '/dev/stdout' (FORMAT SAM, INCLUDE_HEADER false, COMPRESSION gzip);
