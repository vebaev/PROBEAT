#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace="${PROBEAT_ROOT:-$(cd "$script_dir/../.." && pwd)}"
helper_env="${PROBEAT_HELPER_ENV:-$workspace/mamba_envs/snakemake-conda-helper-2026-06-04}"
db_root="${PROBEAT_DATABASES:-$workspace/data/databases}"
downloads="${PROBEAT_DOWNLOADS:-$workspace/downloads/databases}"
DATABASE_ROOT="$db_root"
BAGEL4_COMMIT="${BAGEL4_COMMIT:-dcb961a8cf7487d7affd7788881b49e64b25313b}"
GUTSMASH_COMMIT="${GUTSMASH_COMMIT:-1651919f7e93d09a995ea080a3c029be84ee48c6}"

export PATH="$helper_env/bin:$PATH"
export LANG=C
export LC_ALL=C
export TMPDIR="$workspace/.tmp"
export XDG_CACHE_HOME="$workspace/.cache"
export MPLCONFIGDIR="$workspace/.cache/matplotlib"

mkdir -p "$db_root" "$downloads" "$TMPDIR" "$XDG_CACHE_HOME" "$MPLCONFIGDIR" logs/databases

find_env_bin() {
    local command="$1"
    local required="${2:-}"
    local env_dir
    for env_dir in "$workspace"/snakemake_conda_envs/*_; do
        [[ -x "$env_dir/bin/$command" ]] || continue
        if [[ -n "$required" && ! -x "$env_dir/bin/$required" ]]; then
            continue
        fi
        printf '%s\n' "$env_dir/bin/$command"
        return 0
    done
    echo "Could not find executable '$command' in snakemake_conda_envs" >&2
    exit 2
}

ensure_tool_env_bin() {
    local command="$1"
    local env_name="$2"
    local package_spec="$3"
    local env_dir="$workspace/mamba_envs/$env_name"

    if [[ -x "$env_dir/bin/$command" ]]; then
        printf '%s\n' "$env_dir/bin/$command"
        return 0
    fi

    if [[ -x "$helper_env/bin/mamba" ]]; then
        run_log "${env_name}_env_create" \
            "$helper_env/bin/mamba" create -y -p "$env_dir" \
            -c conda-forge -c bioconda "$package_spec"
    else
        run_log "${env_name}_env_create" \
            "$helper_env/bin/conda" create -y -p "$env_dir" \
            -c conda-forge -c bioconda "$package_spec"
    fi

    if [[ ! -x "$env_dir/bin/$command" ]]; then
        echo "Could not create executable '$command' in $env_dir" >&2
        exit 2
    fi
    printf '%s\n' "$env_dir/bin/$command"
}

find_or_create_env_bin() {
    local command="$1"
    local env_name="$2"
    local package_spec="$3"
    local env_dir
    for env_dir in "$workspace"/snakemake_conda_envs/*_; do
        [[ -x "$env_dir/bin/$command" ]] || continue
        printf '%s\n' "$env_dir/bin/$command"
        return 0
    done
    ensure_tool_env_bin "$command" "$env_name" "$package_spec"
}

ensure_yaml_env_bin() {
    local command="$1"
    local env_name="$2"
    local yaml_path="$3"
    local env_dir="$workspace/mamba_envs/$env_name"

    if [[ -x "$env_dir/bin/$command" ]]; then
        printf '%s\n' "$env_dir/bin/$command"
        return 0
    fi

    if [[ ! -s "$yaml_path" ]]; then
        echo "Missing environment definition: $yaml_path" >&2
        exit 2
    fi

    run_log "${env_name}_env_create" \
        "$helper_env/bin/conda" env create --quiet --no-default-packages \
        --file "$yaml_path" --prefix "$env_dir"

    if [[ ! -x "$env_dir/bin/$command" ]]; then
        echo "Could not create executable '$command' in $env_dir" >&2
        exit 2
    fi
    printf '%s\n' "$env_dir/bin/$command"
}

run_log() {
    local name="$1"
    shift
    echo "==> $name" >&2
    "$@" > "logs/databases/${name}.log" 2>&1
}

tool_env() {
    local cmd="$1"
    local env_dir
    env_dir="$(dirname "$(dirname "$cmd")")"
    env PATH="$(dirname "$cmd"):$PATH" \
        PERL5LIB="$env_dir/lib/perl5/site_perl:$env_dir/lib/perl5/core_perl:${PERL5LIB:-}"
}

install_amrfinderplus() {
    local cmd
    cmd="$(find_env_bin amrfinder_update)"
    mkdir -p "$db_root/amrfinderplus"
    run_log amrfinderplus \
        env PATH="$(dirname "$cmd"):$PATH" "$cmd" \
        --database "$db_root/amrfinderplus" --force_update
}

install_mobsuite() {
    local cmd
    cmd="$(find_env_bin mob_init makeblastdb)"
    mkdir -p "$db_root/mob_suite"
    run_log mob_suite \
        env PATH="$(dirname "$cmd"):$PATH" "$cmd" \
        --database_directory "$db_root/mob_suite"
}

install_busco_bacteria() {
    local cmd
    cmd="$(find_or_create_env_bin busco busco-6.0.0 busco=6.0.0)"
    mkdir -p "$db_root/busco"
    run_log busco_bacteria \
        env PATH="$(dirname "$cmd"):$PATH" "$cmd" \
        --download bacteria_odb10 --download_path "$db_root/busco"
}

install_antismash() {
    local cmd
    cmd="$(find_env_bin download-antismash-databases)"
    mkdir -p "$db_root/antismash"
    run_log antismash \
        env PATH="$(dirname "$cmd"):$PATH" "$cmd" \
        --database-dir "$db_root/antismash"
}

install_bakta_light() {
    local cmd
    cmd="$(find_env_bin bakta_db)"
    local bakta_downloads="$downloads/bakta"
    local extracted="$bakta_downloads/db-light"
    mkdir -p "$db_root/bakta" "$bakta_downloads"

    if [[ -s "$db_root/bakta/version.json" ]]; then
        echo "Bakta light database already installed at $db_root/bakta"
        return 0
    fi

    if [[ ! -s "$extracted/version.json" ]]; then
        run_log bakta_light_download \
            env PATH="$(dirname "$cmd"):$PATH" "$cmd" \
            download --type light --output "$bakta_downloads"
    fi

    if [[ -s "$extracted/version.json" ]]; then
        run_log bakta_light_install cp -a "$extracted/." "$db_root/bakta/"
        return 0
    fi

    local db_file
    db_file="$(find "$bakta_downloads" -type f \( -name '*light*.tar.gz' -o -name '*light*.tar.xz' -o -name '*light*.tar' \) -printf '%T@ %p\n' | sort -nr | awk 'NR == 1 {print $2}')"
    if [[ -z "${db_file:-}" ]]; then
        echo "Could not find downloaded Bakta light database tarball in $bakta_downloads" >&2
        exit 2
    fi

    run_log bakta_light_install \
        env PATH="$(dirname "$cmd"):$PATH" "$cmd" \
        install --db-file "$db_file" --output "$db_root/bakta"
}

install_bakta_full() {
    local cmd
    cmd="$(find_env_bin bakta_db)"
    local bakta_downloads="$downloads/bakta"
    local extracted="$bakta_downloads/db"
    local target_dir="$db_root/bakta_full"
    mkdir -p "$target_dir" "$bakta_downloads"

    if [[ -s "$target_dir/version.json" ]]; then
        echo "Bakta full database already installed at $target_dir"
        return 0
    fi

    if [[ ! -s "$extracted/version.json" ]]; then
        run_log bakta_full_download \
            env PATH="$(dirname "$cmd"):$PATH" "$cmd" \
            download --type full --minor 0 --output "$bakta_downloads"
    fi

    if [[ -s "$extracted/version.json" ]]; then
        run_log bakta_full_install cp -a "$extracted/." "$target_dir/"
        return 0
    fi

    local db_file
    db_file="$(find "$bakta_downloads" -type f \( -name '*full*.tar.gz' -o -name '*full*.tar.xz' -o -name '*full*.tar' -o -name 'db.tar.gz' -o -name 'db.tar.xz' -o -name 'db.tar' \) -printf '%T@ %p\n' | sort -nr | awk 'NR == 1 {print $2}')"
    if [[ -z "${db_file:-}" ]]; then
        echo "Could not find downloaded Bakta full database tarball in $bakta_downloads" >&2
        exit 2
    fi

    run_log bakta_full_install \
        env PATH="$(dirname "$cmd"):$PATH" "$cmd" \
        install --db-file "$db_file" --output "$target_dir"
}

install_gtdbtk() {
    local cmd
    cmd="$(find_env_bin gtdbtk)"
    local gtdb_downloads="$downloads/gtdbtk"
    local archive="$gtdb_downloads/gtdbtk_r226_data.tar.gz"
    local url="https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz"

    mkdir -p "$db_root/gtdbtk" "$gtdb_downloads"

    if [[ ! -s "$archive" ]]; then
        run_log gtdbtk_download curl -L --fail --continue-at - --output "$archive" "$url"
    fi

    if [[ ! -s "$db_root/gtdbtk/metadata/metadata.txt" ]]; then
        run_log gtdbtk_extract tar -xzf "$archive" -C "$db_root/gtdbtk" --strip-components 1
    fi

    run_log gtdbtk_check env GTDBTK_DATA_PATH="$db_root/gtdbtk" PATH="$(dirname "$cmd"):$PATH" "$cmd" check_install
}

install_mash_refseq() {
    local mash_dir="$db_root/mash"
    local archive="$downloads/mash/RefSeqSketchesDefaults.msh.gz"
    local db="$mash_dir/RefSeqSketchesDefaults.msh"
    local host="gembox.cbcb.umd.edu"
    local ip="128.8.120.25"
    local url="https://$host/mash/RefSeqSketchesDefaults.msh.gz"

    mkdir -p "$mash_dir" "$downloads/mash"

    if [[ -s "$db" ]]; then
        echo "Mash RefSeq database already installed at $db"
        return 0
    fi

    if [[ ! -s "$archive" ]]; then
        run_log mash_refseq_download curl -L --fail --continue-at - --resolve "$host:443:$ip" --output "$archive" "$url"
    fi

    run_log mash_refseq_extract bash -c "gzip -dc '$archive' > '$db'"
}

install_sourmash_gtdb_rs220() {
    local cmd
    cmd="$(find_or_create_env_bin sourmash sourmash sourmash)"
    local sourmash_dir="$db_root/sourmash"
    local base_url="https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs220"
    mkdir -p "$sourmash_dir"

    if [[ ! -s "$sourmash_dir/gtdb-reps-rs220-k31.dna.zip" ]]; then
        run_log sourmash_gtdb_reps_rs220 \
            curl -L --fail --continue-at - \
            --output "$sourmash_dir/gtdb-reps-rs220-k31.dna.zip" \
            "$base_url/gtdb-reps-rs220-k31.dna.zip"
    fi
    if [[ ! -s "$sourmash_dir/gtdb-rs220.lineages.csv" ]]; then
        run_log sourmash_gtdb_lineages_rs220 \
            curl -L --fail --continue-at - \
            --output "$sourmash_dir/gtdb-rs220.lineages.csv" \
            "$base_url/gtdb-rs220.lineages.csv"
    fi

    run_log sourmash_gtdb_check "$cmd" sig summarize "$sourmash_dir/gtdb-reps-rs220-k31.dna.zip"
}

install_ncbi_refseq_species() {
    local ncbi_dir="$db_root/ncbi_refseq_species"
    local source="$ncbi_dir/assembly_summary_refseq_historical.txt"
    local representatives="$ncbi_dir/refseq_historical_species_representatives.tsv"
    mkdir -p "$ncbi_dir"

    run_log ncbi_refseq_historical_download \
        curl -L --fail --continue-at - \
        --output "$source" \
        https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq_historical.txt
    run_log ncbi_refseq_species_filter \
        "$helper_env/bin/python" workflow/scripts/filter_ncbi_species_representatives.py \
        "$source" "$representatives"
}

install_go_basic() {
    local go_dir="$db_root/go"
    local target="$go_dir/go-basic.obo"
    mkdir -p "$go_dir"
    if [[ ! -s "$target" ]]; then
        run_log go_basic_download \
            curl -L --fail --continue-at - \
            --output "$target" \
            https://current.geneontology.org/ontology/go-basic.obo
    fi
}

install_abricate_db() {
    local db_name="$1"
    local cmd
    local env_dir
    cmd="$(find_or_create_env_bin abricate abricate-1.0.1 abricate=1.0.1)"
    env_dir="$(dirname "$(dirname "$cmd")")"
    mkdir -p "$db_root/abricate/$db_name"

    if [[ ! -s "$db_root/abricate/$db_name/sequences" && -s "$env_dir/db/$db_name/sequences" ]]; then
        run_log "abricate_${db_name}_copy" cp -a "$env_dir/db/$db_name/." "$db_root/abricate/$db_name/"
    fi

    run_log "abricate_${db_name}" \
        env PATH="$(dirname "$cmd"):$PATH" \
        PERL5LIB="$env_dir/lib/perl5/site_perl:$env_dir/lib/perl5/core_perl:${PERL5LIB:-}" \
        "$cmd" \
        --datadir "$db_root/abricate" --setupdb
    if [[ ! -s "$db_root/abricate/$db_name/sequences" ]]; then
        echo "Could not find ABRicate database '$db_name' in $db_root/abricate" >&2
        exit 2
    fi
}

install_plasmidfinder() {
    install_abricate_db plasmidfinder
}

install_vfdb() {
    install_abricate_db vfdb
}

install_kofamscan() {
    local cmd
    cmd="$(find_or_create_env_bin exec_annotation kofamscan-1.3.0 kofamscan=1.3.0)"
    local kofam_dir="$db_root/kofam"
    local kofam_downloads="$downloads/kofam"
    local profiles_archive="$kofam_downloads/profiles.tar.gz"
    local ko_list_archive="$kofam_downloads/ko_list.gz"
    mkdir -p "$kofam_dir" "$kofam_downloads"

    if [[ ! -s "$profiles_archive" ]]; then
        run_log kofamscan_profiles_download curl -L --fail --continue-at - --output "$profiles_archive" \
            https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
    fi
    if [[ ! -s "$ko_list_archive" ]]; then
        run_log kofamscan_ko_list_download curl -L --fail --continue-at - --output "$ko_list_archive" \
            https://www.genome.jp/ftp/db/kofam/ko_list.gz
    fi

    if [[ ! -d "$kofam_dir/profiles" ]]; then
        run_log kofamscan_profiles_extract tar -xzf "$profiles_archive" -C "$kofam_dir"
    fi
    if [[ ! -s "$kofam_dir/ko_list" ]]; then
        run_log kofamscan_ko_list_extract bash -c "gzip -dc '$ko_list_archive' > '$kofam_dir/ko_list'"
    fi
}

install_keggdecoder() {
    local env_dir="$workspace/mamba_envs/keggdecoder-1.3"
    if [[ ! -x "$env_dir/bin/KEGG-decoder" ]]; then
        if [[ -x "$env_dir/bin/python" ]]; then
            run_log keggdecoder_env_dependencies \
                "$helper_env/bin/mamba" install -y -p "$env_dir" -c conda-forge \
                pandas matplotlib seaborn
        else
            run_log keggdecoder_env_create \
                "$helper_env/bin/mamba" create -y -p "$env_dir" -c conda-forge \
                python=3.10 pip pandas matplotlib seaborn
        fi
        run_log keggdecoder_install \
            "$env_dir/bin/pip" install --no-deps "KEGGDecoder==1.3"
    fi
    "$env_dir/bin/KEGG-decoder" --help >/dev/null
}

install_rasusa() {
    ensure_tool_env_bin rasusa rasusa-2.1.0 rasusa=2.1.0 >/dev/null
}

install_porechop() {
    ensure_tool_env_bin porechop porechop-0.2.4 porechop=0.2.4 >/dev/null
}

install_flye() {
    ensure_tool_env_bin flye flye-2.9.6 flye=2.9.6 >/dev/null
}

install_medaka() {
    local cmd
    cmd="$(ensure_tool_env_bin medaka_consensus medaka medaka)"
    run_log medaka_model_r1041_e82_400bps_hac_v4.1.0 \
        env PATH="$(dirname "$cmd"):$PATH" \
        "$workspace/mamba_envs/medaka/bin/medaka" tools download_models \
        --models r1041_e82_400bps_hac_v4.1.0
}

install_ont_tools() {
    install_porechop
    install_rasusa
    install_flye
    install_medaka
}

install_bagel4() {
    local bagel_dir="$db_root/bagel4"
    local env_dir="$workspace/mamba_envs/bagel4"

    if [[ ! -d "$bagel_dir/.git" ]]; then
        rm -rf "$bagel_dir"
        run_log bagel4_clone git clone https://github.com/annejong/BAGEL4.git "$bagel_dir"
    fi
    run_log bagel4_checkout git -C "$bagel_dir" checkout --detach "$BAGEL4_COMMIT"

    if [[ ! -x "$env_dir/bin/perl" ]]; then
        run_log bagel4_env_create \
            "$helper_env/bin/mamba" create -y -p "$env_dir" -c conda-forge -c bioconda \
            perl perl-json perl-bioperl perl-dbi perl-dbd-sqlite hmmer blast blast-legacy glimmer
    elif [[ ! -f "$env_dir/lib/perl5/site_perl/DBI.pm" && ! -f "$env_dir/lib/perl5/vendor_perl/DBI.pm" ]] || [[ ! -x "$env_dir/bin/blastp" ]]; then
        run_log bagel4_env_dependencies \
            "$helper_env/bin/mamba" install -y -p "$env_dir" -c conda-forge -c bioconda \
            perl-dbi perl-dbd-sqlite blast
    fi

    run_log bagel4_paths \
        bash -c "
            find '$bagel_dir' -type f \\( -name '*.pl' -o -name '*.conf' \\) -print0 |
                xargs -0 sed -i \
                    -e 's#/data/bagel4#$bagel_dir#g' \
                    -e 's#/data/molgentools/lib#$bagel_dir/lib#g'
            sed -i \
                -e 's#^glimmerpath *=.*#glimmerpath = $env_dir/bin#' \
                -e 's#^hmmsearch *=.*#hmmsearch = $env_dir/bin#' \
                '$bagel_dir/bagel4.conf'
            find '$bagel_dir' -maxdepth 1 -type f -name '*.pl' -exec chmod +x {} +
        "

    test -x "$bagel_dir/bagel4_wrapper.pl"
    test -x "$env_dir/bin/hmmsearch"
    test -x "$env_dir/bin/blastall"
}

install_dbcan() {
    local db_root="${DATABASE_ROOT}/dbcan_v5"
    local cmd="$workspace/mamba_envs/dbcan_v5/bin/run_dbcan"
    if [[ ! -x "$cmd" ]]; then
        cmd="$(find_or_create_env_bin run_dbcan dbcan_v5 dbcan=5.2.9)"
    fi
    echo "Downloading dbCAN v5 database..."
    mkdir -p "$db_root"
    run_log dbcan_database env PATH="$(dirname "$cmd"):$PATH" "$cmd" database --db_dir "$db_root" --aws_s3
    echo "dbCAN database installed to $db_root"
}

install_genomad() {
    local db_root="${DATABASE_ROOT}/genomad"
    local cmd
    cmd="$(find_or_create_env_bin genomad genomad genomad)"
    mkdir -p "$db_root"
    echo "Downloading geNomad database..."
    run_log genomad \
        env PATH="$(dirname "$cmd"):$PATH" "$cmd" download-database "$db_root"
    echo "geNomad database installed to $db_root"
}

install_genomespot_models() {
    local db_root="${DATABASE_ROOT}/genomespot"
    if [ -d "$db_root" ] && [ -f "$db_root/instructions.json" ]; then
        echo "GenomeSPOT models already installed at $db_root"
        return 0
    fi
    mkdir -p "$db_root"
    echo "Cloning GenomeSPOT models..."
    local tmpdir=$(mktemp -d)
    git clone --depth 1 https://github.com/cultivarium/GenomeSPOT.git "$tmpdir" 2>/dev/null
    cp -r "$tmpdir/models/"* "$db_root/" 2>/dev/null
    rm -rf "$tmpdir"
    echo "GenomeSPOT models installed to $db_root"
}

install_gutsmash_db() {
    local db_root="${DATABASE_ROOT}/gutsmash"
    local src_dir=""
    local py
    local stage="$downloads/gutsmash/antismash_databases"
    local latest_pfam

    echo "Preparing gutSMASH source..."
    for candidate in "${PROBEAT_GUTSMASH_SRC:-}" "$workspace/workflow/gutsmash_src" "/opt/probeat/gutsmash_src"; do
        if [[ -n "$candidate" && -s "$candidate/run_gutsmash.py" && -s "$candidate/antismash/__main__.py" ]]; then
            src_dir="$candidate"
            break
        fi
    done
    if [[ -z "$src_dir" ]]; then
        src_dir="$workspace/workflow/gutsmash_src"
        if [[ ! -d "$src_dir" ]]; then
            git clone -b gutsmash https://github.com/victoriapascal/gutsmash.git "$src_dir"
        fi
        git -C "$src_dir" checkout --detach "$GUTSMASH_COMMIT"
    fi

    py="$(ensure_yaml_env_bin python gutsmash workflow/envs/gutsmash.yaml)"
    mkdir -p "$db_root" "$stage"

    echo "Downloading and preparing gutSMASH Pfam/resfam/ClusterBlast databases..."
    run_log gutsmash_large_databases \
        env PYTHONPATH="$src_dir:${PYTHONPATH:-}" PATH="$(dirname "$py"):$PATH" \
        "$py" "$src_dir/antismash/download_databases.py" --database-dir "$stage"

    echo "Normalizing gutSMASH database layout..."
    rm -rf "$db_root/.staging"
    mkdir -p "$db_root/.staging"

    latest_pfam="$(find "$stage/pfam" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' 2>/dev/null | sort -V | tail -n 1)"
    if [[ -z "$latest_pfam" || ! -s "$stage/pfam/$latest_pfam/Pfam-A.hmm" ]]; then
        echo "Missing downloaded gutSMASH Pfam database under $stage/pfam" >&2
        exit 2
    fi
    mkdir -p "$db_root/.staging/pfam"
    cp -a "$stage/pfam/$latest_pfam/Pfam-A.hmm"* "$db_root/.staging/pfam/"

    if [[ ! -s "$stage/resfam/Resfams.hmm" ]]; then
        echo "Missing downloaded gutSMASH Resfams database under $stage/resfam" >&2
        exit 2
    fi
    mkdir -p "$db_root/.staging/resfam/resfam"
    cp -a "$stage/resfam/Resfams.hmm"* "$db_root/.staging/resfam/resfam/"

    if [[ ! -s "$stage/clusterblast/clusters.txt" || ! -s "$stage/clusterblast/proteins.fasta" ]]; then
        echo "Missing downloaded gutSMASH ClusterBlast database under $stage/clusterblast" >&2
        exit 2
    fi
    mkdir -p "$db_root/.staging/ClusterBlast"
    cp -a "$stage/clusterblast/clusters.txt" "$stage/clusterblast/proteins.fasta" "$db_root/.staging/ClusterBlast/"
    if [[ -s "$stage/clusterblast/proteins.dmnd" ]]; then
        cp -a "$stage/clusterblast/proteins.dmnd" "$db_root/.staging/ClusterBlast/"
    else
        run_log gutsmash_clusterblast_diamond \
            env PATH="$(dirname "$py"):$PATH" diamond makedb \
            --in "$db_root/.staging/ClusterBlast/proteins.fasta" \
            --db "$db_root/.staging/ClusterBlast/proteins"
    fi

    if [[ ! -s "$src_dir/antismash/modules/clusterblast/data/known/clusters.txt" || ! -s "$src_dir/antismash/modules/clusterblast/data/known/proteins.fasta" ]]; then
        echo "Missing bundled gutSMASH KnownClusterBlast data in $src_dir" >&2
        exit 2
    fi
    mkdir -p "$db_root/.staging/KnownClusterBlast"
    cp -a "$src_dir/antismash/modules/clusterblast/data/known/clusters.txt" \
          "$src_dir/antismash/modules/clusterblast/data/known/proteins.fasta" \
          "$db_root/.staging/KnownClusterBlast/"
    if [[ -s "$src_dir/antismash/modules/clusterblast/data/known/proteins.dmnd" ]]; then
        cp -a "$src_dir/antismash/modules/clusterblast/data/known/proteins.dmnd" "$db_root/.staging/KnownClusterBlast/"
    else
        run_log gutsmash_knownclusterblast_diamond \
            env PATH="$(dirname "$py"):$PATH" diamond makedb \
            --in "$db_root/.staging/KnownClusterBlast/proteins.fasta" \
            --db "$db_root/.staging/KnownClusterBlast/proteins"
    fi

    if [[ ! -s "$src_dir/antismash/detection/genefunctions/gut_data/smcogs.hmm" ]]; then
        echo "Missing bundled gutSMASH smcogs data in $src_dir" >&2
        exit 2
    fi
    mkdir -p "$db_root/.staging/smcogs"
    cp -a "$src_dir/antismash/detection/genefunctions/gut_data/smcogs.hmm"* "$db_root/.staging/smcogs/"
    smcogs_needs_press=0
    for suffix in h3f h3i h3m h3p; do
        if [[ ! -s "$db_root/.staging/smcogs/smcogs.hmm.$suffix" || "$db_root/.staging/smcogs/smcogs.hmm" -nt "$db_root/.staging/smcogs/smcogs.hmm.$suffix" ]]; then
            smcogs_needs_press=1
        fi
    done
    if [[ "$smcogs_needs_press" -eq 1 ]]; then
        run_log gutsmash_smcogs_hmmpress \
            env PATH="$(dirname "$py"):$PATH" hmmpress -f "$db_root/.staging/smcogs/smcogs.hmm"
    fi

    rm -rf "$db_root/pfam" "$db_root/resfam" "$db_root/ClusterBlast" "$db_root/KnownClusterBlast" "$db_root/smcogs"
    mv "$db_root/.staging/pfam" "$db_root/pfam"
    mv "$db_root/.staging/resfam" "$db_root/resfam"
    mv "$db_root/.staging/ClusterBlast" "$db_root/ClusterBlast"
    mv "$db_root/.staging/KnownClusterBlast" "$db_root/KnownClusterBlast"
    mv "$db_root/.staging/smcogs" "$db_root/smcogs"
    rmdir "$db_root/.staging"

    test -s "$db_root/pfam/Pfam-A.hmm"
    test -s "$db_root/resfam/resfam/Resfams.hmm"
    test -s "$db_root/ClusterBlast/clusters.txt"
    test -s "$db_root/ClusterBlast/proteins.fasta"
    test -s "$db_root/ClusterBlast/proteins.dmnd"
    test -s "$db_root/KnownClusterBlast/clusters.txt"
    test -s "$db_root/KnownClusterBlast/proteins.fasta"
    test -s "$db_root/KnownClusterBlast/proteins.dmnd"
    test -s "$db_root/smcogs/smcogs.hmm"
    test -s "$db_root/smcogs/smcogs.hmm.h3f"
    test -s "$db_root/smcogs/smcogs.hmm.h3i"
    test -s "$db_root/smcogs/smcogs.hmm.h3m"
    test -s "$db_root/smcogs/smcogs.hmm.h3p"
    echo "gutSMASH databases installed to $db_root"
}

install_eggnog() {
    local db_root="${DATABASE_ROOT}/eggnog"
    local cmd
    cmd="$(find_or_create_env_bin download_eggnog_data.py eggnog eggnog-mapper=2.1.15)"
    mkdir -p "$db_root"
    echo "Downloading eggNOG v7 database..."
    run_log eggnog \
        env PATH="$(dirname "$cmd"):$PATH" "$cmd" --data_dir "$db_root"
    echo "eggNOG database installed to $db_root"
}

usage() {
    cat <<'USAGE'
Usage: workflow/scripts/install_databases.sh TARGET [TARGET...]

Targets:
  amrfinderplus
  busco-bacteria
  mobsuite
  antismash
  bakta-light
  bakta-full
  mash-refseq    Mash RefSeq sketch database
  sourmash-gtdb-rs220
                 sourmash GTDB RS220 representative k=31 database
  ncbi-refseq-species
                 one historical RefSeq representative assembly per species
  go-basic       Gene Ontology basic OBO file
  plasmidfinder  PlasmidFinder through ABRicate
  vfdb           VFDB through ABRicate
  abricate-probiotic
                 PlasmidFinder and VFDB through ABRicate
  kofamscan      KOfam profiles and ko_list
  keggdecoder    KEGGDecoder tool; uses KofamScan output, no extra database
  rasusa         Rasusa read subsampling tool; no database
  porechop       Porechop ONT adapter trimming tool
  flye           Flye long-read assembler
  medaka         Medaka ONT consensus polishing tool
  ont-tools      Porechop, Rasusa, Flye, and Medaka
  bagel4         BAGEL4 bacteriocin/RiPP mining tool and bundled databases
  dbcan          dbCAN CAZyme database without CGC/PUL data
  genomad        geNomad plasmid and virus identification database
  eggnog         eggNOG v7 functional annotation database
  genomespot     GenomeSPOT growth condition prediction models
  gutsmash       gutSMASH metabolic gene cluster database
  probiotic      plasmidfinder, vfdb, kofamscan, dbcan
  gtdbtk        GTDB-Tk R226 full package; very large
  small          amrfinderplus, busco-bacteria, mobsuite
  all-light      small, mash-refseq, bakta-light, antismash

GTDB-Tk is intentionally not included in all-light because the reference
package is very large and should be downloaded as an explicit storage-planned
step.
USAGE
}

if [[ $# -eq 0 ]]; then
    usage
    exit 2
fi

for target in "$@"; do
    case "$target" in
        amrfinderplus) install_amrfinderplus ;;
        busco-bacteria) install_busco_bacteria ;;
        mobsuite) install_mobsuite ;;
        antismash) install_antismash ;;
        bakta-light) install_bakta_light ;;
        bakta-full) install_bakta_full ;;
        mash-refseq) install_mash_refseq ;;
        sourmash-gtdb-rs220) install_sourmash_gtdb_rs220 ;;
        ncbi-refseq-species) install_ncbi_refseq_species ;;
        go-basic) install_go_basic ;;
        plasmidfinder) install_plasmidfinder ;;
        vfdb) install_vfdb ;;
        abricate-probiotic)
            install_plasmidfinder
            install_vfdb
            ;;
        kofamscan) install_kofamscan ;;
        keggdecoder) install_keggdecoder ;;
        rasusa) install_rasusa ;;
        porechop) install_porechop ;;
        flye) install_flye ;;
        medaka) install_medaka ;;
        ont-tools) install_ont_tools ;;
        bagel4) install_bagel4 ;;
        dbcan) install_dbcan ;;
        genomad) install_genomad ;;
        eggnog) install_eggnog ;;
        genomespot) install_genomespot_models ;;
        gutsmash) install_gutsmash_db ;;
        probiotic)
            install_plasmidfinder
            install_vfdb
            install_kofamscan
            install_dbcan
            ;;
        gtdbtk) install_gtdbtk ;;
        small)
            install_amrfinderplus
            install_busco_bacteria
            install_mobsuite
            ;;
        all-light)
            install_amrfinderplus
            install_busco_bacteria
            install_mobsuite
            install_mash_refseq
            install_bakta_light
            install_antismash
            ;;
        -h|--help) usage ;;
        *)
            echo "Unknown target: $target" >&2
            usage >&2
            exit 2
            ;;
    esac
done
