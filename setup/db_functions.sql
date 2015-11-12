-- Function: refreshview()

-- DROP FUNCTION refreshview();

CREATE OR REPLACE FUNCTION refreshview()
  RETURNS void AS
$BODY$
begin

  execute format('DROP VIEW metadata_view');
  execute format('CREATE OR REPLACE VIEW metadata_view AS SELECT *
FROM  (SELECT g.id, gs.external_id_prefix||''_''||g.id_at_source as genome, g.name as organism_name, g.description, u.username as username
      FROM genomes as g
      LEFT JOIN genome_sources as gs ON gs.id=g.genome_source_id
      LEFT JOIN users as u ON u.id=g.owner_id) as temp
LEFT JOIN metadata_nucleotide as m1 USING (id)
LEFT JOIN metadata_genes as m2 USING (id)
LEFT JOIN metadata_taxonomy as m3 USING (id)
LEFT JOIN metadata_ncbi as m4 USING (id)');

  execute format('ALTER TABLE metadata_view
  OWNER TO gtdb');
end
$BODY$
  LANGUAGE plpgsql VOLATILE SECURITY DEFINER
  COST 100;
ALTER FUNCTION refreshview()
  OWNER TO gtdb;



-- Function: upsert(regclass, character varying, character varying, text[], text[])

-- DROP FUNCTION upsert(regclass, character varying, character varying, text[], text[]);

CREATE OR REPLACE FUNCTION upsert(
    metatable regclass,
    field character varying,
    typemeta character varying,
    _arr1_id_source text[],
    _arr2_meta text[])
  RETURNS void AS
$BODY$

BEGIN
    EXECUTE format('CREATE TEMPORARY TABLE newvals(id integer,id_source text, newdata %s);',typemeta);
    EXECUTE format('INSERT INTO newvals(id_source, newdata) SELECT unnest($1), unnest(cast($2 as %s[]))',typemeta) USING _arr1_id_source,_arr2_meta;

    EXECUTE format('UPDATE newvals
	SET id = g.id
	FROM genomes g
	WHERE g.id_at_source = newvals.id_source;',metatable ,field);

    EXECUTE format('LOCK TABLE %s IN EXCLUSIVE MODE',metatable);

    EXECUTE format('UPDATE %I as meta
	SET %I = newvals.newdata
	FROM newvals
	WHERE newvals.id = meta.id;',metatable ,field);

    EXECUTE format('INSERT INTO %I (id,%I)
	SELECT newvals.id, newvals.newdata
	FROM newvals
	LEFT OUTER JOIN %I ON (%I.id = newvals.id)
	WHERE %I.id IS NULL;',metatable ,field, metatable,metatable,metatable,metatable);
    DROP TABLE newvals;
    EXCEPTION WHEN foreign_key_violation THEN
            RAISE EXCEPTION 'One foreign key is not present in the genome table';
	WHEN invalid_text_representation THEN
	    RAISE EXCEPTION 'Please Verify Value of the TSV file. The Data on the right column can not be cast to the selected datatype';
	WHEN not_null_violation THEN
	    RAISE EXCEPTION 'Genome Ids present in the TSV file is not in GTDB Database';
	WHEN feature_not_supported THEN
	    RAISE EXCEPTION 'Impossible to run the command';

    COMMIT;

END
$BODY$
  LANGUAGE plpgsql VOLATILE
  COST 100;
ALTER FUNCTION upsert(regclass, character varying, character varying, text[], text[])
  OWNER TO gtdb;



-- Function: upsert_aligned_markers(integer, integer, boolean, text, boolean)

-- DROP FUNCTION upsert_aligned_markers(integer, integer, boolean, text, boolean);

CREATE OR REPLACE FUNCTION upsert_aligned_markers(
    g_id integer,
    m_id integer,
    is_dna boolean,
    seq text,
    multi_hit boolean)
  RETURNS void AS
$BODY$
BEGIN
    LOOP
        -- first try to update the key
        UPDATE aligned_markers
            SET sequence = seq,
            multiple_hits = multi_hit
            WHERE genome_id = g_id
            AND marker_id = m_id
            AND dna = is_dna;
        IF found THEN
            RETURN;
        END IF;
        -- not there, so try to insert the key
        -- if someone else inserts the same key concurrently,
        -- we could get a unique-key failure
        BEGIN
            INSERT INTO aligned_markers(genome_id, marker_id, dna, sequence, multiple_hits)
                VALUES (g_id, m_id, is_dna, seq, multi_hit);
            RETURN;
        EXCEPTION WHEN unique_violation THEN
            -- Do nothing, and loop to try the UPDATE again.
        END;
    END LOOP;
END;
$BODY$
  LANGUAGE plpgsql VOLATILE
  COST 100;
ALTER FUNCTION upsert_aligned_markers(integer, integer, boolean, text, boolean)
  OWNER TO gtdb;
