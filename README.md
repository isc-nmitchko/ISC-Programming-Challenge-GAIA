```
       .          *               .        *           .            .
   *        .            ____      _      ___      _          .         *
        .              / ___|    / \    |_ _|    / \      *         .
    *         .       | |  _    / _ \    | |    / _ \         .          *
          .           | |_| |  / ___ \   | |   / ___ \    .        *
      *        .       \____| /_/   \_\ |___| /_/   \_\        .       *
   .        *      .          .            *          .            .    *
        D R 3   E P O C H   P H O T O M E T R Y   .   V A R I A B I L I T Y
  _____________________________________________________________________________
 |  Find every source whose BP or RP flux changed by more than 100% over the   |
 |  Gaia DR3 observation window. InterSystems IRIS + Embedded Python.          |
 |_____________________________________________________________________________|
```

# Gaia DR3 Variable-Source Finder

A submission to the **1st InterSystems Programming Challenge**
([contest #47](https://openexchange.intersystems.com/contest/47)).

The application reads the Gaia DR3 epoch-photometry archive, and for every
`source_id` it computes how much each band's flux varied across all valid
observations:

```
                      max_flux - min_flux
   percentage_change = ------------------- x 100        (computed per band)
                            min_flux

   result = max( percentage_change_BP , percentage_change_RP )
```

Invalid samples (`NaN`, null, empty) are ignored. Every source whose
`percentage_change` exceeds **100%** is written to the output, one record per
line:

```
   source_id , bp_min_flux , bp_max_flux , rp_min_flux , rp_max_flux , percentage_change
```

On the official 20-file benchmark this scans **75,068** sources and reports
**57,099** variable ones.


```
   .-=========================================================================-.
   ||  T W O   E N G I N E S ,   O N E   A N S W E R   ( b y t e - i d e n t ) ||
   '-=========================================================================-'
```

| Engine                       | Entry point                    | Optimized for      | Warm time |
|------------------------------|--------------------------------|--------------------|----------:|
| Polars (default)             | `src/gaia.py`                  | Code Golf, clarity |  ~1.2 s   |
| libdeflate + C scan (opt-in) | `src/gaia_fast.py` + `.c`      | Benchmarking       |  ~1.0 s   |

Both engines produce **identical output to the byte**. The default is a compact,
readable polars pipeline; the opt-in engine drops to a hand-written C kernel for
the fastest honest run we could measure.


## Installation & Running

Prerequisites: [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
and [Docker Desktop](https://www.docker.com/products/docker-desktop).

```bash
git clone <this-repo>
cd intersystems-challenge1-docker-template
docker compose up --build -d
```

The build installs polars into IRIS Embedded Python, pre-compiles the optional C
kernel, and compiles `RunScript.mac`. Run the challenge exactly as the judges do:

```bash
docker compose exec iris iris session iris -U USER
```
```objectscript
USER> do ^RunScript
Matched sources: 57099
Elapsed time: 1.5 seconds
```

The answer is written to `data/out/result.csv`:

```
10655814178816,23.783341359526982,157841.99293824387,35.634526749900566,1566.8893953334514,663566.1794160244
10892037246720,203.7002183789311,666.4385312062207,332.5918636164385,661.7290497218023,227.16633124392905
...
```


## How It Works

The input files are gzip-compressed **ECSV** (a CSV body behind a ~365-line `#`
YAML header). Each row is one source; `bp_flux` and `rp_flux` are quoted arrays
such as `"[1820.8,2013.8,NaN,...]"`. Only three of the ~47 columns matter:
`source_id`, `bp_flux`, `rp_flux`.

### Engine 1 - Polars (default, `src/gaia.py`)

```
   *.csv.gz --(polars Rust gzip)--> read 3 cols --> per-file min/max -.
   *.csv.gz --(polars Rust gzip)--> read 3 cols --> per-file min/max --+--> concat
   *.csv.gz --(polars Rust gzip)--> read 3 cols --> per-file min/max -'      |
                  (16 worker threads, largest file first)                    v
                                          percentage_change -> filter>100 -> CSV
```

Each file is read on its own thread, with polars decompressing gzip natively in
Rust. `NaN` is stripped at the string level so the float cast yields null (which
`list.min`/`list.max` already skip), and each worker reduces its file to five
scalar columns *before* the cross-file concat. Roughly a dozen meaningful lines.

### Engine 2 - libdeflate + C kernel (opt-in, `src/gaia_fast.py`, `src/gaiascan.c`)

```
   .--------------------- gaiascan.so  (OpenMP, GIL released) ---------------------.
   |  *.gz --> libdeflate_gzip_decompress --> single-pass two-column min/max scan  |
   '-------------------------------------------------------------------------------'
                       |  struct-of-arrays  (no per-element Python)
                       v
        polars:  NaN -> null  ->  percentage_change  ->  filter > 100  ->  write_csv
```

For each file, in parallel across all cores with the GIL released, the kernel:

1. Decompresses with **libdeflate** - the fastest gzip decoder available - sizing
   the output buffer from the gzip trailer's `ISIZE` field (no reallocation).
2. Scans only the three needed columns in a **single forward pass**: skip the `#`
   header, count commas (respecting quoted fields), and for the two flux arrays
   walk the values with `strtod`, tracking four running min/max doubles. No
   DataFrame, no per-value allocation, no storage of individual fluxes.
3. Returns results as struct-of-arrays, which Python copies into polars in bulk.
   The percentage-change / filter / CSV-write tail is shared with Engine 1, so the
   output is identical.

It is engineered to never fail: the `.so` is pre-built at image build, lazily
recompiled at first run if the dev volume mount shadows it, and falls back to the
pure-polars engine if `gcc`/`libdeflate` are unavailable.

### Selecting the fast engine

`RunScript.mac` imports the default engine. To benchmark the C kernel, change one
line:

```objectscript
; Set gaia = ##class(%SYS.Python).Import("gaia")        ; polars (default)
  Set gaia = ##class(%SYS.Python).Import("gaia_fast")   ; libdeflate + C
```


## The Optimization Story

This solution was built by measuring iterations. Every change had to be both
**faster** and **identical** to the previous answer, verified each round
against a fixed oracle.

```
   10.87s  ##################################################  serial polars
    2.76s  #############                                       parallel reads (16 threads)
    1.69s  ########                                            polars-native Rust gzip
    1.53s  #######                                             removed per-row list.eval
    1.20s  ######                                              largest-first + per-file reduce
    1.02s  #####                                               C kernel + libdeflate (the floor)
```

- **Round 1 - escape Python's overhead.** The first profile showed 93% of the
  time in Python's `gzip` module not data processing with polars. Parallelizing the reads, then
  letting polars decompress in Rust, then deleting a per-row sub-expression took
  10.87 s to 1.53 s.
- **Round 2 - schedule the gzip deflate.** A single gzip stream cannot be split, so
  the 28 MB file bounds the wall clock; we start it first and reduce each file to
  scalars before concatenating. ~1.2 s. We also *rejected with evidence* the
  drop-in decoders `isal` (1.66 s) and `cramjam` (2.45 s) - both lost to
  polars-native, because routing bytes back through Python costs more than any
  faster decoder saves.
- **Round 3 - find the processing floow.** The profiling revealed  the job is
  **decompression-bound**. We measured decompress-only in C at **1.03 s**, equal
  to the full kernel time - the scan and marshalling are under 0.01 s combined.
  That makes ~1.0 s a hard floor: the time to push 1.54 GB through the fastest
  decoder on this hardware. `OMP_NUM_THREADS=16` (the core count) was optimal, and
  `-march=native` gave no gain (the work lives inside precompiled libdeflate), so
  we kept a portable `-O3` build with no `SIGILL` risk on the judge's CPU.

```
   ,---------------------------------------------------------------------------.
   |  The job is decompression-bound. ~1.0s is the floor, and the measurement  |
   |  proves it: nothing reads 1.54 GB of gzip faster without changing the     |
   |  problem (e.g. precomputing a cache, which we deliberately did not do).    |
   '---------------------------------------------------------------------------'
```


## Correctness

The subtle contract honored by every iteration: `NaN`, null, and empty samples
are invalid; a band with no valid values yields a null min/max that
`max_horizontal` ignores, so a source appears only if a real band swung past
100%. After every change we verified:

```
   wc -l result.csv                         ->  57099
   awk -F, 'NF!=6 || $6<=100' result.csv    ->  (no rows)
   diff <(sort oracle) <(sort result.csv)   ->  0 lines   (byte-identical)
```


## Project Layout

```
   src/
     gaia.py         default engine - pure polars, code-golf friendly
     gaia_fast.py    opt-in engine  - libdeflate + C kernel, polars fallback
     gaiascan.c      the native decompress + single-pass scan kernel
     RunScript.mac   benchmark entry point invoked by  do ^RunScript
   data/
     in/             the 20 EpochPhotometry_*.csv.gz benchmark files
     out/result.csv  the answer - 57,099 variable sources
   Dockerfile        IRIS + Embedded Python; installs polars, builds gaiascan.so
   docker-compose.yml
```


## Tech Stack

- **InterSystems IRIS Community + Embedded Python** - the platform; `%SYS.Python`
  bridges ObjectScript and Python directly.
- **[Polars](https://pola.rs/)** - multithreaded, Rust-powered DataFrames; the
  default engine and the shared output tail of both.
- **[libdeflate](https://github.com/ebiggers/libdeflate) + C / OpenMP** - the
  native fast path.
- **Docker** - one `docker compose up --build` reproduces everything.

```
   "We are made of star-stuff."  -  and, measured carefully, of about one second.
```

## Data & Attribution

Input data: ESA Gaia DR3 epoch photometry. Column definitions:
<https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_photometry/ssec_dm_epoch_photometry.html>.
How to cite Gaia: <https://gea.esac.esa.int/archive/documentation/credits.html>.
