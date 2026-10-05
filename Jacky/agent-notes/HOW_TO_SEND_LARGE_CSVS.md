# How to get Runjia’s large CSVs into this agent

We **do** need these three (≈67 MB total). Do **not** put them on GitHub.

| file | size | put at |
|---|---:|---|
| `guide_index.csv` | 28 MB | `Runjia/data/interim/guide_index.csv` |
| `a_flank.csv` | 36 MB | `Runjia/data/interim/features/a_flank.csv` |
| `g_transcription.csv` | 2.8 MB | `Runjia/data/interim/features/g_transcription.csv` |

## Best options (pick one)

### 1. Google Drive / Dropbox / WeTransfer / Box (easiest)

1. Zip the three files (keep the folder structure if possible).
2. Upload to Drive/Dropbox/WeTransfer.
3. Share a **direct download link** with anyone who has the link.
4. Paste the URL in this chat.

I will `curl`/`wget` it into `Runjia/data/interim/`.

**Drive tip:** a normal share link often needs converting to a direct URL, or use “Anyone with the link” + paste here and I will handle it. WeTransfer direct links usually work as-is.

### 2. Cursor chat attachment

If the chat accepts ~70 MB: attach one zip named e.g. `runjia_interim_pack.zip` in the next message. Cursor saves attachments under uploads; I will unpack them.

If the upload fails on size, use option 1.

### 3. Do **not** use

- GitHub commits / PRs (bloated history; you already hit limits)
- Email to me (I cannot receive email)
- Rebuilding `g_transcription` yourself unless you already have RegulonDB + PRECISE-1K — Runjia is right that transfer is cheaper

## Still needed besides these three (public downloads — you or I can fetch)

From Runjia’s note, these are **not** his to redistribute; we download them:

1. `ecoli_feature_matrix.csv` (~189 MB) — Noshay NAR S2  
2. `NC_000913.2.fasta` + `genes_NC_000913.2.tsv` — preferably from Pasteur `badSeed_public` so coords match  
3. crisprHAL `data/*.csv` + `NC_013716.1.fasta` — for transfer bootstrap  
4. Guo supplementary tables — for orientation / ΔrecA  

Once the three CSVs + the E. coli matrix + genome are here, the Priority A workplan runs can start.

## Already received from this upload

- `within_gene_variance_share` = **0.8451** (most label variance is **within**-gene — good news for Thread 4)
- `base_meta.json`, `qct_lookup_tables.json`, `v_column_mapping.json`
- His READ_ME + MANIFEST (saved under `Jacky/agent-notes/from-runjia/`)

Also noted: he asked that `workplan.py` not live inside `Runjia/sgrna/`. Moved to `Jacky/agent-notes/scripts/workplan.py`.
