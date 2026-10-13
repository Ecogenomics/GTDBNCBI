# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Internal (non-public) software for the Genome Taxonomy Database: a PostgreSQL-backed store of NCBI and user genomes, their metadata, marker-gene alignments and taxonomy, plus tooling to build multiple sequence alignments / trees from marker sets. The code was recently migrated from Python 2 to Python 3 (see `Release.txt`), and that migration is still in progress in the working tree.

## CURRENT_TASK.md

Keep `CURRENT_TASK.md` in the repo root up to date with what we are working on right now, so work can resume if the server crashes or the connection drops mid-session. Read it at the start of every session. Update it whenever the task, plan or state changes, such as after each meaningful step, commit, push or PR, not just at the end. Record:

- the goal
- the branch and any open PR
- what is done
- what is in progress, including uncommitted changes
- the next steps
- any open questions or decisions waiting on the user

When a task is finished, say so in the file. The file is gitignored, and this repo is public, so never commit it.

## Setup and commands

- `src/gtdb/Config.py` is gitignored and must exist; create it from `src/gtdb/ConfigTemplate.py`. It holds DB credentials, genome root directories, `DB_SERVERS` (database name → host) and `LATEST_DB`. Never commit it.
- Install: `pip install .` (or `python setup.py install`). `setup.py` installs `bin/gtdb` plus **every file** under `scripts_user/` and `scripts_dev/` as an executable script, and aborts if two scripts share a basename. Keep script names unique across those trees.
- Run: `gtdb [-r] [-u USER] [-db RELEASE] [-t THREADS] <category> <command> ...` (categories: `users`, `genomes`, `genome_lists`, `markers`, `marker_sets`, `metadata`, `taxonomy`, `tree`, `db_stats`, `power`). `gtdb <category> -h` lists the commands in a category.
- There is no test suite, linter config or CI. To check syntax: `python3 -m py_compile bin/gtdb && python3 -m compileall -q src`. Running any real command needs a live PostgreSQL connection.
- `src/gtdb/VERSION` holds both the version and the changelogs. Its first line must be `version=X.Y.Z`: `setup.py` (regex `^version=`) and `gtdb/__init__.py` (first line) both read it to get `__version__`. Below that line are the changelogs: "GTDB software versions", "GTDB database versions" and "Taxonomy versions". `bin/gtdb --version` reports the software version plus `Config.LATEST_DB`. `Release.txt` holds old 2015 notes and is not maintained.
- **Versioning:** follow semantic versioning (MAJOR.MINOR.PATCH). Every change to the software must bump the `version=` line in `src/gtdb/VERSION` and add a dated entry at the top of the "GTDB software versions" changelog in the same file:
  - PATCH: bug fixes.
  - MINOR: new backwards-compatible functionality, such as new commands or flags.
  - MAJOR: incompatible changes, such as CLI changes that break existing usage or code that no longer works with the previous DB schema.

  DB schema changes also get an entry under "GTDB database versions".
- Python dependencies: `biolib < 0.1.0`, `psycopg2`, numpy, etc. External binaries are checked lazily with `biolib.external.execute.check_dependencies` in the `__init__` of the class that needs them (e.g. `Prodigal` → prodigal, `TigrfamSearch`/`PfamSearch` → hmmsearch, tree creation → FastTree/FastTreeMP). Follow this pattern; don't check them globally.

## Architecture

**Layering:** `bin/gtdb` (argparse CLI) → small handler function per command → method on the `GenomeDatabase` facade (`src/gtdb/GenomeDatabase.py`) → `*Manager` classes.

- **`bin/gtdb`** builds the parsers, connects with `db.conn.MakePostgresConnection(args.release)`, calls `db.Login(...)`, then dispatches through `args.func(db, args)`. A falsy return value means failure, and the CLI then dumps `db.GetErrors()` / `db.GetWarnings()`.
- **`GenomeDatabase`** owns the psycopg2 connection (`GenomeDatabaseConnection`), the current `User`, and lists of errors and warnings. Each public method opens a cursor, creates the managers it needs as `XManager(cur, self.currentUser[, threads])`, and owns the transaction: `self.conn.commit()` on success; on `GenomeDatabaseError`, `self.ReportError(...)` and `return False`; on any other exception, `rollback()` and re-raise. Managers don't commit; they raise `GenomeDatabaseError` for expected failures.
- **Managers:** `GenomeManager` (adding, deleting and pulling genomes, which runs Prodigal, TIGRFAM and Pfam with `multiprocessing` and moves files into the GTDB directory layout), `GenomeListManager` and `MarkerSetManager` (named sets with owners and private flags), `MarkerManager`, `AlignedMarkerManager` (hmmalign per genome, cached in the `aligned_markers` table through the `upsert_aligned_markers` SQL function), `TreeManager` (genome filtering by quality and taxa, MSA column trimming, and writing ARB/metadata/FASTA outputs), `MetadataManager`, `GenomeRepresentativeManager`, `PowerUserManager` (admin reports and checks), `UserManager`.
- **Auth:** login uses the Linux username (`pwd`) matched against the `users` table. `-r` logs in as root (the user must have root rights), and `-u USER` impersonates a user (this needs root).
- **Release locking:** `-db` picks any database in `Config.DB_SERVERS`. Anything other than `Config.LATEST_DB` sets `GenomeDatabase.dbLock`, and while it is set, add and delete operations refuse to run. `Config.DB_UPDATE = True` also blocks add and delete during maintenance.
- **Schema:** `setup/gtdb_create.sql` (tables, plus the stored functions the managers call). `setup/new_postgres_db.sql` bootstraps a fresh database.
- **Output files:** commands that write a single `--output` file open it with `Tools.openOutputFile()`, which gzip-compresses it when the path ends in `.gz`. New commands should do the same, and their help text should mention it.
- **Constants:** thresholds, reference DB paths and per-genome file suffixes (e.g. `_protein.faa.gz`, `_tigrfam_lite_tophit.tsv.gz`) live in `DefaultValues.py` / `ConfigMetadata.py`. Genome files on disk are gzipped.
- **Genome IDs:** external IDs are prefixed `GB_`/`RS_` (GenBank/RefSeq NCBI genomes) or `U_` (user genomes). `TreeManager.canonical_gid` handles this canonicalisation.

**`scripts_dev/`** holds standalone maintenance scripts, mostly for loading metadata into the DB at each release, scraping LPSN, DSMZ and StrainInfo, and handling NCBI FTP/rsync updates. Scripts that use the DB follow `scripts_dev/template_db_script.py`: they instantiate `GenomeDatabase`, log in, and run raw SQL through `db.conn.cursor()`. Several subdirectories have their own gitignored `database_configuration/Config*Update.py`. `legacy_code/`, `deprecated/` and `deprecated_scripts/` are unmaintained. `docs/how_to/` has release procedures (NCBI rsync, taxonomy update).

## Python 3 migration gotchas

- Use `str(e)` for exception messages. Python 3 exceptions have no `.message`, and `GenomeDatabaseError` doesn't define one.
- Most scripts in `scripts_dev/` are still Python 2 only. About 62 of the 86 outside `legacy_code/` and `deprecated*/` don't parse under Python 3, mostly because of `print` statements. Port a script before relying on it. `template_db_script.py` has already been ported.
- Check for other py2 idioms in any code you touch: `iteritems`, `raw_input`, binary vs. text file modes (`'wt'` for text output), and integer `/`.
- Stale `*.pyc` files from Python 2 sit next to sources in `src/gtdb/`. Ignore them.
