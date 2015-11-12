--
-- PostgreSQL database dump
--

-- Dumped from database version 9.1.19
-- Dumped by pg_dump version 9.4.0
-- Started on 2015-11-12 11:08:27

SET statement_timeout = 0;
SET lock_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SET check_function_bodies = false;
SET client_min_messages = warning;

--
-- TOC entry 194 (class 3079 OID 11677)
-- Name: plpgsql; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;


--
-- TOC entry 2094 (class 0 OID 0)
-- Dependencies: 194
-- Name: EXTENSION plpgsql; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION plpgsql IS 'PL/pgSQL procedural language';


SET search_path = public, pg_catalog;

--
-- TOC entry 208 (class 1255 OID 95184)
-- Name: refreshview(); Type: FUNCTION; Schema: public; Owner: gtdb
--

CREATE FUNCTION refreshview() RETURNS void
    LANGUAGE plpgsql SECURITY DEFINER
    AS $$
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
LEFT JOIN metadata_ncbi as m4 USING (id)
LEFT JOIN lpsn_genera as lpsn1 ON m3.gtdb_genus=lpsn1.lpsn_genus
LEFT JOIN lpsn_species as lpsn2 ON m3.gtdb_species=lpsn2.lpsn_species
LEFT JOIN lpsn_strains as lpsn3 ON m4.ncbi_organism_name=lpsn3.lpsn_strain');

  execute format('ALTER TABLE metadata_view
  OWNER TO gtdb');
end
$$;


ALTER FUNCTION public.refreshview() OWNER TO gtdb;

--
-- TOC entry 206 (class 1255 OID 95185)
-- Name: upsert(regclass, character varying, character varying, text[], text[]); Type: FUNCTION; Schema: public; Owner: gtdb
--

CREATE FUNCTION upsert(metatable regclass, field character varying, typemeta character varying, _arr1_id_source text[], _arr2_meta text[]) RETURNS void
    LANGUAGE plpgsql
    AS $_$

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
$_$;


ALTER FUNCTION public.upsert(metatable regclass, field character varying, typemeta character varying, _arr1_id_source text[], _arr2_meta text[]) OWNER TO gtdb;

--
-- TOC entry 207 (class 1255 OID 95186)
-- Name: upsert_aligned_markers(integer, integer, boolean, text, boolean); Type: FUNCTION; Schema: public; Owner: gtdb
--

CREATE FUNCTION upsert_aligned_markers(g_id integer, m_id integer, is_dna boolean, seq text, multi_hit boolean) RETURNS void
    LANGUAGE plpgsql
    AS $$
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
$$;


ALTER FUNCTION public.upsert_aligned_markers(g_id integer, m_id integer, is_dna boolean, seq text, multi_hit boolean) OWNER TO gtdb;

SET default_tablespace = '';

SET default_with_oids = false;

--
-- TOC entry 167 (class 1259 OID 95187)
-- Name: aligned_markers; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE aligned_markers (
    genome_id integer NOT NULL,
    marker_id integer NOT NULL,
    sequence text,
    multiple_hits boolean NOT NULL,
    evalue text,
    bitscore text
);


ALTER TABLE aligned_markers OWNER TO gtdb;

--
-- TOC entry 168 (class 1259 OID 95193)
-- Name: genome_list_contents; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE genome_list_contents (
    list_id integer NOT NULL,
    genome_id integer NOT NULL
);


ALTER TABLE genome_list_contents OWNER TO gtdb;

--
-- TOC entry 169 (class 1259 OID 95196)
-- Name: genome_lists; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE genome_lists (
    id integer NOT NULL,
    name text NOT NULL,
    description text,
    owned_by_root boolean DEFAULT false NOT NULL,
    owner_id integer,
    private boolean DEFAULT true,
    display_order integer DEFAULT 1000
);


ALTER TABLE genome_lists OWNER TO gtdb;

--
-- TOC entry 2095 (class 0 OID 0)
-- Dependencies: 169
-- Name: COLUMN genome_lists.display_order; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN genome_lists.display_order IS 'Indicates desired order of items when displayed to user.';


--
-- TOC entry 170 (class 1259 OID 95205)
-- Name: genome_lists_id_seq; Type: SEQUENCE; Schema: public; Owner: gtdb
--

CREATE SEQUENCE genome_lists_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE genome_lists_id_seq OWNER TO gtdb;

--
-- TOC entry 2096 (class 0 OID 0)
-- Dependencies: 170
-- Name: genome_lists_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: gtdb
--

ALTER SEQUENCE genome_lists_id_seq OWNED BY genome_lists.id;


--
-- TOC entry 171 (class 1259 OID 95207)
-- Name: genome_sources; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE genome_sources (
    id integer NOT NULL,
    name text NOT NULL,
    external_id_prefix text NOT NULL,
    last_auto_id integer DEFAULT 0 NOT NULL,
    user_editable boolean DEFAULT false NOT NULL
);


ALTER TABLE genome_sources OWNER TO gtdb;

--
-- TOC entry 172 (class 1259 OID 95215)
-- Name: genome_sources_id_seq; Type: SEQUENCE; Schema: public; Owner: gtdb
--

CREATE SEQUENCE genome_sources_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE genome_sources_id_seq OWNER TO gtdb;

--
-- TOC entry 2097 (class 0 OID 0)
-- Dependencies: 172
-- Name: genome_sources_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: gtdb
--

ALTER SEQUENCE genome_sources_id_seq OWNED BY genome_sources.id;


--
-- TOC entry 173 (class 1259 OID 95217)
-- Name: genomes; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE genomes (
    id integer NOT NULL,
    name text NOT NULL,
    description text,
    owned_by_root boolean DEFAULT false NOT NULL,
    owner_id integer,
    fasta_file_location text NOT NULL,
    fasta_file_sha256 text NOT NULL,
    genome_source_id integer NOT NULL,
    id_at_source text NOT NULL,
    date_added timestamp without time zone NOT NULL,
    has_changed boolean DEFAULT true NOT NULL,
    last_update date,
    genes_file_location text,
    genes_file_sha256 text
);


ALTER TABLE genomes OWNER TO gtdb;

--
-- TOC entry 174 (class 1259 OID 95225)
-- Name: genomes_id_seq; Type: SEQUENCE; Schema: public; Owner: gtdb
--

CREATE SEQUENCE genomes_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE genomes_id_seq OWNER TO gtdb;

--
-- TOC entry 2098 (class 0 OID 0)
-- Dependencies: 174
-- Name: genomes_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: gtdb
--

ALTER SEQUENCE genomes_id_seq OWNED BY genomes.id;


--
-- TOC entry 190 (class 1259 OID 95471)
-- Name: lpsn_genera; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE lpsn_genera (
    lpsn_genus text NOT NULL,
    lpsn_type_genus text,
    lpsn_genus_authority text
);


ALTER TABLE lpsn_genera OWNER TO gtdb;

--
-- TOC entry 2099 (class 0 OID 0)
-- Dependencies: 190
-- Name: TABLE lpsn_genera; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON TABLE lpsn_genera IS 'Genus information from LPSN.';


--
-- TOC entry 2100 (class 0 OID 0)
-- Dependencies: 190
-- Name: COLUMN lpsn_genera.lpsn_genus; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN lpsn_genera.lpsn_genus IS 'Genus name.';


--
-- TOC entry 2101 (class 0 OID 0)
-- Dependencies: 190
-- Name: COLUMN lpsn_genera.lpsn_type_genus; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN lpsn_genera.lpsn_type_genus IS 'Indicates if genus is a type genus of a family.';


--
-- TOC entry 2102 (class 0 OID 0)
-- Dependencies: 190
-- Name: COLUMN lpsn_genera.lpsn_genus_authority; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN lpsn_genera.lpsn_genus_authority IS 'Reference information provided by LPSN.';


--
-- TOC entry 191 (class 1259 OID 95482)
-- Name: lpsn_species; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE lpsn_species (
    lpsn_species text NOT NULL,
    lpsn_type_species text,
    lpsn_species_authority text
);


ALTER TABLE lpsn_species OWNER TO gtdb;

--
-- TOC entry 2103 (class 0 OID 0)
-- Dependencies: 191
-- Name: TABLE lpsn_species; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON TABLE lpsn_species IS 'Species information from LPSN.';


--
-- TOC entry 2104 (class 0 OID 0)
-- Dependencies: 191
-- Name: COLUMN lpsn_species.lpsn_species; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN lpsn_species.lpsn_species IS 'Species name.';


--
-- TOC entry 2105 (class 0 OID 0)
-- Dependencies: 191
-- Name: COLUMN lpsn_species.lpsn_type_species; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN lpsn_species.lpsn_type_species IS 'Indicates if species is a type species of a genus.';


--
-- TOC entry 2106 (class 0 OID 0)
-- Dependencies: 191
-- Name: COLUMN lpsn_species.lpsn_species_authority; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN lpsn_species.lpsn_species_authority IS 'Reference information provided by LPSN.';


--
-- TOC entry 192 (class 1259 OID 95494)
-- Name: lpsn_strains; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE lpsn_strains (
    lpsn_strain text NOT NULL
);


ALTER TABLE lpsn_strains OWNER TO gtdb;

--
-- TOC entry 2107 (class 0 OID 0)
-- Dependencies: 192
-- Name: TABLE lpsn_strains; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON TABLE lpsn_strains IS 'Strain information from LPSN.';


--
-- TOC entry 2108 (class 0 OID 0)
-- Dependencies: 192
-- Name: COLUMN lpsn_strains.lpsn_strain; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN lpsn_strains.lpsn_strain IS 'Strain name.';


--
-- TOC entry 175 (class 1259 OID 95227)
-- Name: marker_databases; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE marker_databases (
    id integer NOT NULL,
    name text NOT NULL,
    external_id_prefix text NOT NULL,
    last_auto_id integer DEFAULT 0 NOT NULL,
    user_editable boolean DEFAULT false NOT NULL
);


ALTER TABLE marker_databases OWNER TO gtdb;

--
-- TOC entry 176 (class 1259 OID 95235)
-- Name: marker_databases_id_seq; Type: SEQUENCE; Schema: public; Owner: gtdb
--

CREATE SEQUENCE marker_databases_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE marker_databases_id_seq OWNER TO gtdb;

--
-- TOC entry 2109 (class 0 OID 0)
-- Dependencies: 176
-- Name: marker_databases_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: gtdb
--

ALTER SEQUENCE marker_databases_id_seq OWNED BY marker_databases.id;


--
-- TOC entry 177 (class 1259 OID 95237)
-- Name: marker_set_contents; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE marker_set_contents (
    set_id integer NOT NULL,
    marker_id integer NOT NULL
);


ALTER TABLE marker_set_contents OWNER TO gtdb;

--
-- TOC entry 178 (class 1259 OID 95240)
-- Name: marker_sets; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE marker_sets (
    id integer NOT NULL,
    name text NOT NULL,
    description text,
    owned_by_root boolean DEFAULT false NOT NULL,
    owner_id integer,
    private boolean DEFAULT true
);


ALTER TABLE marker_sets OWNER TO gtdb;

--
-- TOC entry 179 (class 1259 OID 95248)
-- Name: marker_sets_id_seq; Type: SEQUENCE; Schema: public; Owner: gtdb
--

CREATE SEQUENCE marker_sets_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE marker_sets_id_seq OWNER TO gtdb;

--
-- TOC entry 2110 (class 0 OID 0)
-- Dependencies: 179
-- Name: marker_sets_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: gtdb
--

ALTER SEQUENCE marker_sets_id_seq OWNED BY marker_sets.id;


--
-- TOC entry 180 (class 1259 OID 95250)
-- Name: markers; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE markers (
    id integer NOT NULL,
    name text NOT NULL,
    description text,
    owned_by_root boolean DEFAULT false NOT NULL,
    owner_id integer,
    marker_file_location text NOT NULL,
    marker_file_sha256 text NOT NULL,
    marker_database_id integer NOT NULL,
    id_in_database text NOT NULL,
    size integer NOT NULL
);


ALTER TABLE markers OWNER TO gtdb;

--
-- TOC entry 181 (class 1259 OID 95257)
-- Name: markers_id_seq; Type: SEQUENCE; Schema: public; Owner: gtdb
--

CREATE SEQUENCE markers_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE markers_id_seq OWNER TO gtdb;

--
-- TOC entry 2111 (class 0 OID 0)
-- Dependencies: 181
-- Name: markers_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: gtdb
--

ALTER SEQUENCE markers_id_seq OWNED BY markers.id;


--
-- TOC entry 182 (class 1259 OID 95259)
-- Name: metadata_genes; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE metadata_genes (
    id integer NOT NULL,
    checkm_completeness double precision,
    checkm_contamination double precision,
    protein_count integer,
    coding_bases integer,
    coding_density double precision,
    ssu_count integer,
    checkm_marker_count integer,
    checkm_marker_lineage text,
    checkm_genome_count integer,
    checkm_marker_set_count integer,
    checkm_strain_heterogeneity double precision
);


ALTER TABLE metadata_genes OWNER TO gtdb;

--
-- TOC entry 2112 (class 0 OID 0)
-- Dependencies: 182
-- Name: COLUMN metadata_genes.id; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_genes.id IS 'id is equivalent to genome_id. 
id is used during the VIEW LEFT JOIN and need to have similar name between the metadata table and the genome table';


--
-- TOC entry 2113 (class 0 OID 0)
-- Dependencies: 182
-- Name: COLUMN metadata_genes.checkm_completeness; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_genes.checkm_completeness IS 'Estimated completeness of genome.';


--
-- TOC entry 2114 (class 0 OID 0)
-- Dependencies: 182
-- Name: COLUMN metadata_genes.checkm_contamination; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_genes.checkm_contamination IS 'Estimated contamination of genome.';


--
-- TOC entry 2115 (class 0 OID 0)
-- Dependencies: 182
-- Name: COLUMN metadata_genes.protein_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_genes.protein_count IS 'Number of protein coding genes.';


--
-- TOC entry 2116 (class 0 OID 0)
-- Dependencies: 182
-- Name: COLUMN metadata_genes.coding_bases; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_genes.coding_bases IS 'Number of coding bases in genome.';


--
-- TOC entry 2117 (class 0 OID 0)
-- Dependencies: 182
-- Name: COLUMN metadata_genes.coding_density; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_genes.coding_density IS 'Percentage of coding bases in genome.';


--
-- TOC entry 2118 (class 0 OID 0)
-- Dependencies: 182
-- Name: COLUMN metadata_genes.ssu_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_genes.ssu_count IS 'Number of 16S rRNA genes identified in genome.';


--
-- TOC entry 2119 (class 0 OID 0)
-- Dependencies: 182
-- Name: COLUMN metadata_genes.checkm_marker_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_genes.checkm_marker_count IS 'Number of marker genes used to estimate genome quality.';


--
-- TOC entry 2120 (class 0 OID 0)
-- Dependencies: 182
-- Name: COLUMN metadata_genes.checkm_marker_lineage; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_genes.checkm_marker_lineage IS 'Lineage or taxonomic group used to determine marker genes';


--
-- TOC entry 2121 (class 0 OID 0)
-- Dependencies: 182
-- Name: COLUMN metadata_genes.checkm_genome_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_genes.checkm_genome_count IS 'Number of genomes used to determine marker genes.';


--
-- TOC entry 2122 (class 0 OID 0)
-- Dependencies: 182
-- Name: COLUMN metadata_genes.checkm_marker_set_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_genes.checkm_marker_set_count IS 'Number of marker sets used to estimate genome quality.';


--
-- TOC entry 2123 (class 0 OID 0)
-- Dependencies: 182
-- Name: COLUMN metadata_genes.checkm_strain_heterogeneity; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_genes.checkm_strain_heterogeneity IS 'Estimated strain heterogeneity of genome.';


--
-- TOC entry 183 (class 1259 OID 95265)
-- Name: metadata_ncbi; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE metadata_ncbi (
    id integer NOT NULL,
    ncbi_biosample text,
    ncbi_total_gap_length integer,
    ncbi_molecule_count integer,
    ncbi_date text,
    ncbi_submitter text,
    ncbi_ncrna_count integer,
    ncbi_scaffold_n50 integer,
    ncbi_assembly_name text,
    ncbi_scaffold_n75 integer,
    ncbi_protein_count integer,
    ncbi_assembly_type text,
    ncbi_rrna_count integer,
    ncbi_genbank_assembly_accession text,
    ncbi_total_length integer,
    ncbi_unspanned_gaps integer,
    ncbi_taxid integer,
    ncbi_trna_count integer,
    ncbi_genome_representation text,
    ncbi_top_level_count integer,
    ncbi_spanned_gaps integer,
    ncbi_translation_table integer,
    ncbi_scaffold_n90 integer,
    ncbi_contig_count integer,
    ncbi_organism_name text,
    ncbi_region_count integer,
    ncbi_contig_n50 integer,
    ncbi_ungapped_length integer,
    ncbi_scaffold_l50 integer,
    ncbi_ssu_count integer,
    ncbi_scaffold_count integer,
    ncbi_assembly_level text,
    ncbi_refseq_assembly_and_genbank_assemblies_identical text,
    ncbi_release_type text
);


ALTER TABLE metadata_ncbi OWNER TO gtdb;

--
-- TOC entry 2124 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.id; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.id IS 'id is equivalent to genome_id. 
id is used during the VIEW LEFT JOIN and need to have similar name between the metadata table and the genome table';


--
-- TOC entry 2125 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_biosample; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_biosample IS 'BioSample identifier.';


--
-- TOC entry 2126 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_total_gap_length; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_total_gap_length IS 'Total length of gaps.';


--
-- TOC entry 2127 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_molecule_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_molecule_count IS 'Number of chromosomes and plasmids in full assembly.';


--
-- TOC entry 2128 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_date; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_date IS '"Date assembly was released in the INSDC (DDBJ, ENA, or GenBank) databases."';


--
-- TOC entry 2129 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_submitter; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_submitter IS '"Individual, institution, or consortium that submitted the assembly."';


--
-- TOC entry 2130 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_ncrna_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_ncrna_count IS 'Number of ncRNA genes identified in genome.';


--
-- TOC entry 2131 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_scaffold_n50; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_scaffold_n50 IS 'Scaffold length at which 50% of total bases in assembly are in scaffolds of that length or greater.';


--
-- TOC entry 2132 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_assembly_name; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_assembly_name IS 'Name of assembly specified by submitter.';


--
-- TOC entry 2133 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_scaffold_n75; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_scaffold_n75 IS 'Scaffold length at which 75% of total bases in assembly are in contigs of that length or greater.';


--
-- TOC entry 2134 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_protein_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_protein_count IS 'Number of protein coding genes.';


--
-- TOC entry 2135 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_assembly_type; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_assembly_type IS 'Ploidy information for genome.';


--
-- TOC entry 2136 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_rrna_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_rrna_count IS 'Number of rRNA genes identified in genome.';


--
-- TOC entry 2137 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_genbank_assembly_accession; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_genbank_assembly_accession IS 'GenBank assembly accesion.';


--
-- TOC entry 2138 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_total_length; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_total_length IS 'Total sequence length including bases and gaps.';


--
-- TOC entry 2139 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_unspanned_gaps; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_unspanned_gaps IS '"Number of unspanned gaps (i.e., gaps between scaffolds)."';


--
-- TOC entry 2140 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_taxid; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_taxid IS 'NCBI taxonomy identifier.';


--
-- TOC entry 2141 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_trna_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_trna_count IS 'Number of tRNA genes identified in genome.';


--
-- TOC entry 2142 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_genome_representation; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_genome_representation IS 'Flag indicating if the goal of the assembly was to represent the complete genome (full) or only part of it (partial).';


--
-- TOC entry 2143 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_top_level_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_top_level_count IS '"Number of chromosomes or plasmids, unplaced/unlocalized scaffolds, alt_loci scaffolds, and patch scaffolds."';


--
-- TOC entry 2144 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_spanned_gaps; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_spanned_gaps IS 'Number of spanned gaps (i.e. gaps within a scaffold).';


--
-- TOC entry 2145 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_translation_table; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_translation_table IS 'Specifies the genetic code used for translation.';


--
-- TOC entry 2146 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_scaffold_n90; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_scaffold_n90 IS 'scaffold length at which 90% of total bases in assembly are in contigs of that length or greater.';


--
-- TOC entry 2147 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_contig_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_contig_count IS 'Number of contigs.';


--
-- TOC entry 2148 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_organism_name; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_organism_name IS 'Name of organism specified by submitter.';


--
-- TOC entry 2149 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_region_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_region_count IS 'Number of genomic regions defined in full assembly.';


--
-- TOC entry 2150 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_contig_n50; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_contig_n50 IS 'Contig length at which 50% of total bases in assembly are in contigs of that length or greater.';


--
-- TOC entry 2151 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_ungapped_length; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_ungapped_length IS 'Total length excluding gaps.';


--
-- TOC entry 2152 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_scaffold_l50; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_scaffold_l50 IS '"Number of scaffolds that are longer than, or equal to, the N50 length."';


--
-- TOC entry 2153 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_ssu_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_ssu_count IS 'Number of 16S rRNA genes identified in genome.';


--
-- TOC entry 2154 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_scaffold_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_scaffold_count IS 'Number of scaffolds.';


--
-- TOC entry 2155 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_assembly_level; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_assembly_level IS '"Highest level of assembly for any object in the assembly (complete genome, scaffold, or contig)."';


--
-- TOC entry 2156 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_refseq_assembly_and_genbank_assemblies_identical; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_refseq_assembly_and_genbank_assemblies_identical IS 'Flag indicating if RefSeq and GenBank assemblies are identical.';


--
-- TOC entry 2157 (class 0 OID 0)
-- Dependencies: 183
-- Name: COLUMN metadata_ncbi.ncbi_release_type; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_ncbi.ncbi_release_type IS '"Specifies if current assembly was a major, minor, or patch release."';


--
-- TOC entry 184 (class 1259 OID 95271)
-- Name: metadata_nucleotide; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE metadata_nucleotide (
    id integer NOT NULL,
    scaffold_count integer,
    gc_count integer,
    longest_scaffold integer,
    gc_percentage double precision,
    total_gap_length integer,
    genome_size integer,
    n50_contigs integer,
    n50_scaffolds integer,
    l50_scaffolds integer,
    contig_count integer,
    ambiguous_bases integer,
    longest_contig integer,
    l50_contigs integer,
    mean_scaffold_length integer,
    mean_contig_length integer
);


ALTER TABLE metadata_nucleotide OWNER TO gtdb;

--
-- TOC entry 2158 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.id; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.id IS 'id is equivalent to genome_id. 
id is used during the VIEW LEFT JOIN and need to have similar name between the metadata table and the genome table';


--
-- TOC entry 2159 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.scaffold_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.scaffold_count IS 'Number of scaffolds in genome.';


--
-- TOC entry 2160 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.gc_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.gc_count IS 'Number of G or C bases in genome.';


--
-- TOC entry 2161 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.longest_scaffold; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.longest_scaffold IS 'Number of bases in longest scaffold.';


--
-- TOC entry 2162 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.gc_percentage; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.gc_percentage IS 'GC content of genome.';


--
-- TOC entry 2163 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.total_gap_length; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.total_gap_length IS 'Number of ambiguous bases comprising gaps in scaffolds.';


--
-- TOC entry 2164 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.genome_size; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.genome_size IS '"Total base pairs in genome including nucleotide bases, ambiguous bases, and gaps."';


--
-- TOC entry 2165 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.n50_contigs; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.n50_contigs IS 'Contig length at which 50% of total bases in assembly are in contigs of that length or greater.';


--
-- TOC entry 2166 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.n50_scaffolds; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.n50_scaffolds IS 'Scaffold length at which 50% of total bases in assembly are in scaffolds of that length or greater.';


--
-- TOC entry 2167 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.l50_scaffolds; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.l50_scaffolds IS '"Number of scaffolds longer than, or equal to, the scaffold N50 length."';


--
-- TOC entry 2168 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.contig_count; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.contig_count IS 'Number of contigs in genome.';


--
-- TOC entry 2169 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.ambiguous_bases; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.ambiguous_bases IS 'Number of ambiguous bases in contigs.';


--
-- TOC entry 2170 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.longest_contig; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.longest_contig IS 'Number of bases in longest contig.';


--
-- TOC entry 2171 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.l50_contigs; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.l50_contigs IS '"Number of contigs longer than, or equal to, the contig N50 length."';


--
-- TOC entry 2172 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.mean_scaffold_length; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.mean_scaffold_length IS 'Mean length of scaffolds in base pairs.';


--
-- TOC entry 2173 (class 0 OID 0)
-- Dependencies: 184
-- Name: COLUMN metadata_nucleotide.mean_contig_length; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_nucleotide.mean_contig_length IS 'Mean length of contigs in base pairs.';


--
-- TOC entry 185 (class 1259 OID 95274)
-- Name: metadata_taxonomy; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE metadata_taxonomy (
    id integer NOT NULL,
    ssu_gg_2013_08_taxonomy text,
    ssu_gg_2013_08_blast_bitscore integer,
    ssu_gg_2013_08_blast_subject_id text,
    ssu_gg_2013_08_blast_perc_identity double precision,
    ssu_gg_2013_08_blast_evalue double precision,
    ssu_gg_2013_08_blast_align_len integer,
    ssu_gg_2013_08_query_id text,
    ssu_gg_2013_08_length integer,
    ssu_silva_199_gg_taxa_blast_bitscore integer,
    ssu_silva_199_gg_taxa_taxonomy text,
    ssu_silva_199_gg_taxa_blast_subject_id text,
    ssu_silva_199_gg_taxa_query_id text,
    ssu_silva_199_gg_taxa_blast_align_len integer,
    ssu_silva_199_gg_taxa_length integer,
    ssu_silva_199_gg_taxa_blast_evalue double precision,
    ssu_silva_199_gg_taxa_blast_perc_identity double precision,
    ncbi_type_strain text,
    ncbi_taxonomy text,
    gtdb_class text,
    gtdb_species text,
    gtdb_taxonomy text,
    gtdb_phylum text,
    gtdb_family text,
    gtdb_domain text,
    gtdb_order text,
    gtdb_genus text
);


ALTER TABLE metadata_taxonomy OWNER TO gtdb;

--
-- TOC entry 2174 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.id; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.id IS 'id is equivalent to genome_id. 
id is used during the VIEW LEFT JOIN and need to have similar name between the metadata table and the genome table';


--
-- TOC entry 2175 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_gg_2013_08_taxonomy; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_gg_2013_08_taxonomy IS 'Greengenes taxonomy string.';


--
-- TOC entry 2176 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_gg_2013_08_blast_bitscore; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_gg_2013_08_blast_bitscore IS 'Bitscore of BLAST hit.';


--
-- TOC entry 2177 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_gg_2013_08_blast_subject_id; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_gg_2013_08_blast_subject_id IS 'Sequence identifier of BLAST hit.';


--
-- TOC entry 2178 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_gg_2013_08_blast_perc_identity; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_gg_2013_08_blast_perc_identity IS 'Percent identity of BLAST hit.';


--
-- TOC entry 2179 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_gg_2013_08_blast_evalue; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_gg_2013_08_blast_evalue IS 'E-value of BLAST hit.';


--
-- TOC entry 2180 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_gg_2013_08_blast_align_len; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_gg_2013_08_blast_align_len IS 'Alignment length of BLAST hit.';


--
-- TOC entry 2181 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_gg_2013_08_query_id; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_gg_2013_08_query_id IS 'Query id.';


--
-- TOC entry 2182 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_gg_2013_08_length; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_gg_2013_08_length IS 'Length of query sequence.';


--
-- TOC entry 2183 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_blast_bitscore; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_blast_bitscore IS 'Bitscore of BLAST hit.';


--
-- TOC entry 2184 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_taxonomy; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_taxonomy IS 'Greengenes taxonomy string mapped to SILVA database.';


--
-- TOC entry 2185 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_blast_subject_id; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_blast_subject_id IS 'Sequence identifier of BLAST hit.';


--
-- TOC entry 2186 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_query_id; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_query_id IS 'Query id.';


--
-- TOC entry 2187 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_blast_align_len; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_blast_align_len IS 'Alignment length of BLAST hit.';


--
-- TOC entry 2188 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_length; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_length IS 'Length of query sequence.';


--
-- TOC entry 2189 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_blast_evalue; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_blast_evalue IS 'E-value of BLAST hit.';


--
-- TOC entry 2190 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_blast_perc_identity; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ssu_silva_199_gg_taxa_blast_perc_identity IS 'Percent identity of BLAST hit.';


--
-- TOC entry 2191 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ncbi_type_strain; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ncbi_type_strain IS 'NCBI type material.';


--
-- TOC entry 2192 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.ncbi_taxonomy; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.ncbi_taxonomy IS 'NCBI taxonomy string.';


--
-- TOC entry 2193 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.gtdb_class; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.gtdb_class IS 'Class in genome tree taxonomy.';


--
-- TOC entry 2194 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.gtdb_species; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.gtdb_species IS 'Species in genome tree taxonomy.';


--
-- TOC entry 2195 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.gtdb_taxonomy; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.gtdb_taxonomy IS 'Genome tree taxonomy string.';


--
-- TOC entry 2196 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.gtdb_phylum; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.gtdb_phylum IS 'Phylum in genome tree taxonomy.';


--
-- TOC entry 2197 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.gtdb_family; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.gtdb_family IS 'Family in genome tree taxonomy.';


--
-- TOC entry 2198 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.gtdb_domain; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.gtdb_domain IS 'Domain in genome tree taxonomy.';


--
-- TOC entry 2199 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.gtdb_order; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.gtdb_order IS 'Order in genome tree taxonomy.';


--
-- TOC entry 2200 (class 0 OID 0)
-- Dependencies: 185
-- Name: COLUMN metadata_taxonomy.gtdb_genus; Type: COMMENT; Schema: public; Owner: gtdb
--

COMMENT ON COLUMN metadata_taxonomy.gtdb_genus IS 'Genus in genome tree taxonomy.';


--
-- TOC entry 186 (class 1259 OID 95280)
-- Name: users; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE users (
    id integer NOT NULL,
    username text NOT NULL,
    role_id integer NOT NULL,
    has_root_login boolean DEFAULT false NOT NULL
);


ALTER TABLE users OWNER TO gtdb;

--
-- TOC entry 193 (class 1259 OID 95669)
-- Name: metadata_view; Type: VIEW; Schema: public; Owner: gtdb
--

CREATE VIEW metadata_view AS
SELECT temp.id, temp.genome, temp.organism_name, temp.description, temp.username, m1.scaffold_count, m1.gc_count, m1.longest_scaffold, m1.gc_percentage, m1.total_gap_length, m1.genome_size, m1.n50_contigs, m1.n50_scaffolds, m1.l50_scaffolds, m1.contig_count, m1.ambiguous_bases, m1.longest_contig, m1.l50_contigs, m1.mean_scaffold_length, m1.mean_contig_length, m2.checkm_completeness, m2.checkm_contamination, m2.protein_count, m2.coding_bases, m2.coding_density, m2.ssu_count, m2.checkm_marker_count, m2.checkm_marker_lineage, m2.checkm_genome_count, m2.checkm_marker_set_count, m2.checkm_strain_heterogeneity, m3.ssu_gg_2013_08_taxonomy, m3.ssu_gg_2013_08_blast_bitscore, m3.ssu_gg_2013_08_blast_subject_id, m3.ssu_gg_2013_08_blast_perc_identity, m3.ssu_gg_2013_08_blast_evalue, m3.ssu_gg_2013_08_blast_align_len, m3.ssu_gg_2013_08_query_id, m3.ssu_gg_2013_08_length, m3.ssu_silva_199_gg_taxa_blast_bitscore, m3.ssu_silva_199_gg_taxa_taxonomy, m3.ssu_silva_199_gg_taxa_blast_subject_id, m3.ssu_silva_199_gg_taxa_query_id, m3.ssu_silva_199_gg_taxa_blast_align_len, m3.ssu_silva_199_gg_taxa_length, m3.ssu_silva_199_gg_taxa_blast_evalue, m3.ssu_silva_199_gg_taxa_blast_perc_identity, m3.ncbi_type_strain, m3.ncbi_taxonomy, m3.gtdb_class, m3.gtdb_species, m3.gtdb_taxonomy, m3.gtdb_phylum, m3.gtdb_family, m3.gtdb_domain, m3.gtdb_order, m3.gtdb_genus, m4.ncbi_biosample, m4.ncbi_total_gap_length, m4.ncbi_molecule_count, m4.ncbi_date, m4.ncbi_submitter, m4.ncbi_ncrna_count, m4.ncbi_scaffold_n50, m4.ncbi_assembly_name, m4.ncbi_scaffold_n75, m4.ncbi_protein_count, m4.ncbi_assembly_type, m4.ncbi_rrna_count, m4.ncbi_genbank_assembly_accession, m4.ncbi_total_length, m4.ncbi_unspanned_gaps, m4.ncbi_taxid, m4.ncbi_trna_count, m4.ncbi_genome_representation, m4.ncbi_top_level_count, m4.ncbi_spanned_gaps, m4.ncbi_translation_table, m4.ncbi_scaffold_n90, m4.ncbi_contig_count, m4.ncbi_organism_name, m4.ncbi_region_count, m4.ncbi_contig_n50, m4.ncbi_ungapped_length, m4.ncbi_scaffold_l50, m4.ncbi_ssu_count, m4.ncbi_scaffold_count, m4.ncbi_assembly_level, m4.ncbi_refseq_assembly_and_genbank_assemblies_identical, m4.ncbi_release_type, lpsn1.lpsn_genus, lpsn1.lpsn_type_genus, lpsn1.lpsn_genus_authority, lpsn2.lpsn_species, lpsn2.lpsn_type_species, lpsn2.lpsn_species_authority, lpsn3.lpsn_strain FROM ((((((((SELECT g.id, ((gs.external_id_prefix || '_'::text) || g.id_at_source) AS genome, g.name AS organism_name, g.description, u.username FROM ((genomes g LEFT JOIN genome_sources gs ON ((gs.id = g.genome_source_id))) LEFT JOIN users u ON ((u.id = g.owner_id)))) temp LEFT JOIN metadata_nucleotide m1 USING (id)) LEFT JOIN metadata_genes m2 USING (id)) LEFT JOIN metadata_taxonomy m3 USING (id)) LEFT JOIN metadata_ncbi m4 USING (id)) LEFT JOIN lpsn_genera lpsn1 ON ((m3.gtdb_genus = lpsn1.lpsn_genus))) LEFT JOIN lpsn_species lpsn2 ON ((m3.gtdb_species = lpsn2.lpsn_species))) LEFT JOIN lpsn_strains lpsn3 ON ((m4.ncbi_organism_name = lpsn3.lpsn_strain)));


ALTER TABLE metadata_view OWNER TO gtdb;

--
-- TOC entry 187 (class 1259 OID 95292)
-- Name: user_roles; Type: TABLE; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE TABLE user_roles (
    id integer NOT NULL,
    name text NOT NULL
);


ALTER TABLE user_roles OWNER TO gtdb;

--
-- TOC entry 188 (class 1259 OID 95298)
-- Name: users_id_seq; Type: SEQUENCE; Schema: public; Owner: gtdb
--

CREATE SEQUENCE users_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE users_id_seq OWNER TO gtdb;

--
-- TOC entry 2201 (class 0 OID 0)
-- Dependencies: 188
-- Name: users_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: gtdb
--

ALTER SEQUENCE users_id_seq OWNED BY users.id;


--
-- TOC entry 189 (class 1259 OID 95300)
-- Name: view_list_meta_columns; Type: VIEW; Schema: public; Owner: gtdb
--

CREATE VIEW view_list_meta_columns AS
SELECT c.table_name AS "table", c.column_name AS field, pgd.description FROM ((pg_statio_all_tables st JOIN pg_description pgd ON ((pgd.objoid = st.relid))) JOIN information_schema.columns c ON ((((((pgd.objsubid = (c.ordinal_position)::integer) AND ((c.table_schema)::name = st.schemaname)) AND ((c.table_name)::name = st.relname)) AND ((c.table_name)::text ~~ 'metadata%'::text)) AND ((c.column_name)::text !~~ 'id'::text))));


ALTER TABLE view_list_meta_columns OWNER TO gtdb;

--
-- TOC entry 1898 (class 2604 OID 95305)
-- Name: id; Type: DEFAULT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY genome_lists ALTER COLUMN id SET DEFAULT nextval('genome_lists_id_seq'::regclass);


--
-- TOC entry 1901 (class 2604 OID 95306)
-- Name: id; Type: DEFAULT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY genome_sources ALTER COLUMN id SET DEFAULT nextval('genome_sources_id_seq'::regclass);


--
-- TOC entry 1904 (class 2604 OID 95307)
-- Name: id; Type: DEFAULT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY genomes ALTER COLUMN id SET DEFAULT nextval('genomes_id_seq'::regclass);


--
-- TOC entry 1907 (class 2604 OID 95308)
-- Name: id; Type: DEFAULT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY marker_databases ALTER COLUMN id SET DEFAULT nextval('marker_databases_id_seq'::regclass);


--
-- TOC entry 1910 (class 2604 OID 95309)
-- Name: id; Type: DEFAULT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY marker_sets ALTER COLUMN id SET DEFAULT nextval('marker_sets_id_seq'::regclass);


--
-- TOC entry 1912 (class 2604 OID 95310)
-- Name: id; Type: DEFAULT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY markers ALTER COLUMN id SET DEFAULT nextval('markers_id_seq'::regclass);


--
-- TOC entry 1914 (class 2604 OID 95311)
-- Name: id; Type: DEFAULT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY users ALTER COLUMN id SET DEFAULT nextval('users_id_seq'::regclass);


--
-- TOC entry 1917 (class 2606 OID 95313)
-- Name: aligned_markers_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY aligned_markers
    ADD CONSTRAINT aligned_markers_pkey PRIMARY KEY (genome_id, marker_id);


--
-- TOC entry 1919 (class 2606 OID 95315)
-- Name: genome_list_contents_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY genome_list_contents
    ADD CONSTRAINT genome_list_contents_pkey PRIMARY KEY (list_id, genome_id);


--
-- TOC entry 1921 (class 2606 OID 95317)
-- Name: genome_lists_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY genome_lists
    ADD CONSTRAINT genome_lists_pkey PRIMARY KEY (id);


--
-- TOC entry 1923 (class 2606 OID 95319)
-- Name: genome_sources_external_id_prefix_key; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY genome_sources
    ADD CONSTRAINT genome_sources_external_id_prefix_key UNIQUE (external_id_prefix);


--
-- TOC entry 1925 (class 2606 OID 95321)
-- Name: genome_sources_name_key; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY genome_sources
    ADD CONSTRAINT genome_sources_name_key UNIQUE (name);


--
-- TOC entry 1927 (class 2606 OID 95323)
-- Name: genome_sources_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY genome_sources
    ADD CONSTRAINT genome_sources_pkey PRIMARY KEY (id);


--
-- TOC entry 1929 (class 2606 OID 95325)
-- Name: genomes_genome_source_id_id_at_source_key; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY genomes
    ADD CONSTRAINT genomes_genome_source_id_id_at_source_key UNIQUE (genome_source_id, id_at_source);


--
-- TOC entry 1931 (class 2606 OID 95327)
-- Name: genomes_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY genomes
    ADD CONSTRAINT genomes_pkey PRIMARY KEY (id);


--
-- TOC entry 1963 (class 2606 OID 95481)
-- Name: lpsn_genera_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY lpsn_genera
    ADD CONSTRAINT lpsn_genera_pkey PRIMARY KEY (lpsn_genus);


--
-- TOC entry 1965 (class 2606 OID 95491)
-- Name: lpsn_species_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY lpsn_species
    ADD CONSTRAINT lpsn_species_pkey PRIMARY KEY (lpsn_species);


--
-- TOC entry 1967 (class 2606 OID 95501)
-- Name: lpsn_strains_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY lpsn_strains
    ADD CONSTRAINT lpsn_strains_pkey PRIMARY KEY (lpsn_strain);


--
-- TOC entry 1933 (class 2606 OID 95329)
-- Name: marker_databases_external_id_prefix_key; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY marker_databases
    ADD CONSTRAINT marker_databases_external_id_prefix_key UNIQUE (external_id_prefix);


--
-- TOC entry 1935 (class 2606 OID 95331)
-- Name: marker_databases_name_key; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY marker_databases
    ADD CONSTRAINT marker_databases_name_key UNIQUE (name);


--
-- TOC entry 1937 (class 2606 OID 95333)
-- Name: marker_databases_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY marker_databases
    ADD CONSTRAINT marker_databases_pkey PRIMARY KEY (id);


--
-- TOC entry 1939 (class 2606 OID 95335)
-- Name: marker_set_contents_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY marker_set_contents
    ADD CONSTRAINT marker_set_contents_pkey PRIMARY KEY (set_id, marker_id);


--
-- TOC entry 1941 (class 2606 OID 95337)
-- Name: marker_sets_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY marker_sets
    ADD CONSTRAINT marker_sets_pkey PRIMARY KEY (id);


--
-- TOC entry 1943 (class 2606 OID 95339)
-- Name: markers_id_in_database_marker_database_id_key; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY markers
    ADD CONSTRAINT markers_id_in_database_marker_database_id_key UNIQUE (id_in_database, marker_database_id);


--
-- TOC entry 1945 (class 2606 OID 95341)
-- Name: markers_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY markers
    ADD CONSTRAINT markers_pkey PRIMARY KEY (id);


--
-- TOC entry 1949 (class 2606 OID 95343)
-- Name: meta_ncbi_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY metadata_ncbi
    ADD CONSTRAINT meta_ncbi_pkey PRIMARY KEY (id);


--
-- TOC entry 1951 (class 2606 OID 95345)
-- Name: meta_nuc_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY metadata_nucleotide
    ADD CONSTRAINT meta_nuc_pkey PRIMARY KEY (id);


--
-- TOC entry 1947 (class 2606 OID 95347)
-- Name: meta_prot_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY metadata_genes
    ADD CONSTRAINT meta_prot_pkey PRIMARY KEY (id);


--
-- TOC entry 1953 (class 2606 OID 95349)
-- Name: meta_tax_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY metadata_taxonomy
    ADD CONSTRAINT meta_tax_pkey PRIMARY KEY (id);


--
-- TOC entry 1959 (class 2606 OID 95351)
-- Name: user_roles_name_key; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY user_roles
    ADD CONSTRAINT user_roles_name_key UNIQUE (name);


--
-- TOC entry 1961 (class 2606 OID 95353)
-- Name: user_roles_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY user_roles
    ADD CONSTRAINT user_roles_pkey PRIMARY KEY (id);


--
-- TOC entry 1955 (class 2606 OID 95355)
-- Name: users_pkey; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY users
    ADD CONSTRAINT users_pkey PRIMARY KEY (id);


--
-- TOC entry 1957 (class 2606 OID 95357)
-- Name: users_username_key; Type: CONSTRAINT; Schema: public; Owner: gtdb; Tablespace: 
--

ALTER TABLE ONLY users
    ADD CONSTRAINT users_username_key UNIQUE (username);


--
-- TOC entry 1915 (class 1259 OID 95358)
-- Name: aligned_markerid_idx; Type: INDEX; Schema: public; Owner: gtdb; Tablespace: 
--

CREATE INDEX aligned_markerid_idx ON aligned_markers USING btree (marker_id);


--
-- TOC entry 1968 (class 2606 OID 95359)
-- Name: aligned_markers_genome_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY aligned_markers
    ADD CONSTRAINT aligned_markers_genome_id_fkey FOREIGN KEY (genome_id) REFERENCES genomes(id) ON UPDATE CASCADE;


--
-- TOC entry 1969 (class 2606 OID 95364)
-- Name: aligned_markers_marker_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY aligned_markers
    ADD CONSTRAINT aligned_markers_marker_id_fkey FOREIGN KEY (marker_id) REFERENCES markers(id) ON UPDATE CASCADE;


--
-- TOC entry 1970 (class 2606 OID 95369)
-- Name: genome_list_contents_genome_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY genome_list_contents
    ADD CONSTRAINT genome_list_contents_genome_id_fkey FOREIGN KEY (genome_id) REFERENCES genomes(id) ON UPDATE CASCADE;


--
-- TOC entry 1971 (class 2606 OID 95374)
-- Name: genome_list_contents_list_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY genome_list_contents
    ADD CONSTRAINT genome_list_contents_list_id_fkey FOREIGN KEY (list_id) REFERENCES genome_lists(id) ON UPDATE CASCADE;


--
-- TOC entry 1972 (class 2606 OID 95379)
-- Name: genome_lists_owner_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY genome_lists
    ADD CONSTRAINT genome_lists_owner_id_fkey FOREIGN KEY (owner_id) REFERENCES users(id) ON UPDATE CASCADE;


--
-- TOC entry 1973 (class 2606 OID 95384)
-- Name: genomes_genome_source_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY genomes
    ADD CONSTRAINT genomes_genome_source_id_fkey FOREIGN KEY (genome_source_id) REFERENCES genome_sources(id) ON UPDATE CASCADE;


--
-- TOC entry 1974 (class 2606 OID 95389)
-- Name: genomes_owner_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY genomes
    ADD CONSTRAINT genomes_owner_id_fkey FOREIGN KEY (owner_id) REFERENCES users(id) ON UPDATE CASCADE;


--
-- TOC entry 1975 (class 2606 OID 95394)
-- Name: marker_set_contents_marker_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY marker_set_contents
    ADD CONSTRAINT marker_set_contents_marker_id_fkey FOREIGN KEY (marker_id) REFERENCES markers(id) ON UPDATE CASCADE;


--
-- TOC entry 1976 (class 2606 OID 95399)
-- Name: marker_set_contents_set_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY marker_set_contents
    ADD CONSTRAINT marker_set_contents_set_id_fkey FOREIGN KEY (set_id) REFERENCES marker_sets(id) ON UPDATE CASCADE;


--
-- TOC entry 1977 (class 2606 OID 95404)
-- Name: marker_sets_owner_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY marker_sets
    ADD CONSTRAINT marker_sets_owner_id_fkey FOREIGN KEY (owner_id) REFERENCES users(id) ON UPDATE CASCADE;


--
-- TOC entry 1978 (class 2606 OID 95409)
-- Name: markers_marker_database_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY markers
    ADD CONSTRAINT markers_marker_database_id_fkey FOREIGN KEY (marker_database_id) REFERENCES marker_databases(id) ON UPDATE CASCADE;


--
-- TOC entry 1980 (class 2606 OID 95414)
-- Name: meta_ncbi_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY metadata_ncbi
    ADD CONSTRAINT meta_ncbi_fkey FOREIGN KEY (id) REFERENCES genomes(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- TOC entry 1981 (class 2606 OID 95419)
-- Name: meta_nuc_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY metadata_nucleotide
    ADD CONSTRAINT meta_nuc_fkey FOREIGN KEY (id) REFERENCES genomes(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- TOC entry 1979 (class 2606 OID 95424)
-- Name: meta_prot_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY metadata_genes
    ADD CONSTRAINT meta_prot_fkey FOREIGN KEY (id) REFERENCES genomes(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- TOC entry 1982 (class 2606 OID 95429)
-- Name: meta_tax_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY metadata_taxonomy
    ADD CONSTRAINT meta_tax_fkey FOREIGN KEY (id) REFERENCES genomes(id) ON UPDATE CASCADE ON DELETE CASCADE;


--
-- TOC entry 1983 (class 2606 OID 95434)
-- Name: users_role_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gtdb
--

ALTER TABLE ONLY users
    ADD CONSTRAINT users_role_id_fkey FOREIGN KEY (role_id) REFERENCES user_roles(id) ON UPDATE CASCADE;


--
-- TOC entry 2093 (class 0 OID 0)
-- Dependencies: 7
-- Name: public; Type: ACL; Schema: -; Owner: postgres
--

REVOKE ALL ON SCHEMA public FROM PUBLIC;
REVOKE ALL ON SCHEMA public FROM postgres;
GRANT ALL ON SCHEMA public TO postgres;
GRANT ALL ON SCHEMA public TO PUBLIC;


--
-- TOC entry 2202 (class 0 OID 0)
-- Dependencies: 189
-- Name: view_list_meta_columns; Type: ACL; Schema: public; Owner: gtdb
--

REVOKE ALL ON TABLE view_list_meta_columns FROM PUBLIC;
REVOKE ALL ON TABLE view_list_meta_columns FROM gtdb;
GRANT ALL ON TABLE view_list_meta_columns TO gtdb;
GRANT SELECT ON TABLE view_list_meta_columns TO PUBLIC;


-- Completed on 2015-11-12 11:08:40

--
-- PostgreSQL database dump complete
--

