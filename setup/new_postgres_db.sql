CREATE table user_roles (
    id integer PRIMARY KEY,
    name text UNIQUE NOT NULL
);

INSERT INTO user_roles VALUES (0, 'admin');
INSERT INTO user_roles VALUES (1, 'user');

CREATE TABLE users (
    id serial PRIMARY KEY,
    username text UNIQUE NOT NULL,
    role_id integer NOT NULL,
    has_root_login boolean NOT NULL default false,
    FOREIGN KEY (role_id) REFERENCES user_roles(id) ON UPDATE CASCADE
);

INSERT INTO users (username, role_id, has_root_login) VALUES ('uqaskars', 0, true);

CREATE TABLE genome_sources (
    id serial PRIMARY KEY,
    name text UNIQUE NOT NULL,
    external_id_prefix text UNIQUE NOT NULL,
    last_auto_id integer NOT NULL DEFAULT 0,
    user_accessible boolean NOT NULL DEFAULT false
);

INSERT INTO genome_sources (name, external_id_prefix, user_accessible) VALUES ('user', 'U', true);
INSERT INTO genome_sources (name, external_id_prefix) VALUES ('IMG', 'IMG');
INSERT INTO genome_sources (name, external_id_prefix) VALUES ('NCBI', 'NCBI');

CREATE TABLE genomes (
    id serial PRIMARY KEY,
    name text NOT NULL,
    description text,
    owned_by_root boolean NOT NULL DEFAULT false,
    owner_id integer,
    fasta_file_location text NOT NULL,
    fasta_file_sha256 text NOT NULL,
    genome_source_id integer NOT NULL,
    id_at_source text NOT NULL,
    date_added timestamp NOT NULL,
    checkm_completeness float NOT NULL,
    checkm_contamination float NOT NULL,
    UNIQUE (genome_source_id, id_at_source),
    FOREIGN KEY (genome_source_id) REFERENCES genome_sources(id) ON UPDATE CASCADE,
    FOREIGN KEY (owner_id) REFERENCES users(id) ON UPDATE CASCADE
);

CREATE TABLE genome_lists (
    id serial PRIMARY KEY,
    name text NOT NULL,
    description text,
    owned_by_root boolean NOT NULL DEFAULT false,
    owner_id integer,
    private bool NOT NULL DEFAULT TRUE,
    FOREIGN KEY (owner_id) REFERENCES users(id) ON UPDATE CASCADE
);

CREATE TABLE genome_list_contents (
    list_id integer NOT NULL,
    genome_id integer NOT NULL,
    PRIMARY KEY (list_id, genome_id),
    FOREIGN KEY (genome_id) REFERENCES genomes(id) ON UPDATE CASCADE,
    FOREIGN KEY (list_id) REFERENCES genome_lists(id) ON UPDATE CASCADE
);

CREATE TABLE marker_databases (
    id serial PRIMARY KEY,
    name text UNIQUE NOT NULL,
    external_id_prefix text UNIQUE NOT NULL,
    user_accessible boolean NOT NULL DEFAULT false
);


INSERT INTO marker_databases (name, external_id_prefix, user_accessible) VALUES ('user', 'M', true);
INSERT INTO marker_databases (name, external_id_prefix) VALUES ('PFAM', 'PFAM');
INSERT INTO marker_databases (name, external_id_prefix) VALUES ('TIGRFAM', 'TIGR');

CREATE TABLE markers (
    id serial PRIMARY KEY,
    name text NOT NULL,
    owned_by_root boolean NOT NULL DEFAULT false,
    owner_id integer,
    marker_file_location text NOT NULL,
    marker_file_sha256 text NOT NULL,
    database_id integer NOT NULL,
    database_specific_id text NOT NULL,
    size integer NOT NULL,
    UNIQUE (database_specific_id, database_id),
    FOREIGN KEY (database_id) REFERENCES marker_databases(id) ON UPDATE CASCADE
);

CREATE TABLE aligned_markers (
    genome_id integer NOT NULL,
    marker_id integer NOT NULL,
    dna boolean NOT NULL,
    sequence text,
    PRIMARY KEY (genome_id, marker_id, dna),
    FOREIGN KEY (marker_id) REFERENCES markers(id) ON UPDATE CASCADE,
    FOREIGN KEY (genome_id) REFERENCES genomes(id) ON UPDATE CASCADE
);

CREATE TABLE marker_sets (
    id serial PRIMARY KEY,
    name text NOT NULL,
    description text,
    owned_by_root boolean NOT NULL DEFAULT false,
    owner_id integer,
    private bool NOT NULL DEFAULT TRUE,
    FOREIGN KEY (owner_id) REFERENCES users(id) ON UPDATE CASCADE
);

CREATE TABLE marker_set_contents (
    set_id integer NOT NULL,
    marker_id integer NOT NULL,
    PRIMARY KEY (set_id, marker_id),
    FOREIGN KEY (marker_id) REFERENCES markers(id) ON UPDATE CASCADE,
    FOREIGN KEY (set_id) REFERENCES marker_sets(id) ON UPDATE CASCADE
);



