"""Gaia DR3 processor -- OPT-IN fast path (libdeflate + C scan kernel).

This is the Benchmarking-nomination variant. It is decompression-bound at the
libdeflate throughput floor (~1.0 s for the 20-file set), ~11% faster warm than
the pure-polars default in gaia.py.

How it works: src/gaiascan.c decompresses each .csv.gz with libdeflate and scans
ONLY the two flux columns in a single pass, in parallel across files (OpenMP, GIL
released). Per-source min/max come back as struct-of-arrays and are marshalled
into polars with bulk copies (zero per-element Python). Output is byte-identical
to gaia.py's polars path.

To use it from RunScript.mac, import "gaia_fast" instead of "gaia":
    Set gaia = ##class(%SYS.Python).Import("gaia_fast")
    Set count = gaia.run()

If gaiascan.so is missing it is lazily compiled on first run; if gcc/libdeflate
are unavailable it transparently falls back to the pure-polars implementation.
"""
import array, ctypes, glob, os, subprocess
import polars as pl

_HERE = os.path.dirname(__file__)
_SO = os.path.join(_HERE, "gaiascan.so")
_SRC = os.path.join(_HERE, "gaiascan.c")
_LIBDEFLATE = "/lib/x86_64-linux-gnu/libdeflate.so.0"


def _ensure_so():
    """Build gaiascan.so on first use if absent (the dev volume mount can shadow
    the image-baked .so). Best-effort: returns True if the .so is usable."""
    if os.path.exists(_SO):
        return True
    if not (os.path.exists(_SRC) and os.path.exists(_LIBDEFLATE)):
        return False
    try:
        subprocess.run(
            ["gcc", "-O3", "-fopenmp", "-shared", "-fPIC", _SRC,
             "-o", _SO, _LIBDEFLATE, "-lm"],
            check=True, capture_output=True)
        return os.path.exists(_SO)
    except Exception:
        return False


def _finalize(df, outpath):
    """Shared tail: NaN (empty array -> no valid value) becomes null so it is
    excluded, compute percentage_change, filter > 100, write the 6-column CSV."""
    out = (df.with_columns(pl.col("bp_min", "bp_max", "rp_min", "rp_max").fill_nan(None))
           .with_columns(percentage_change=pl.max_horizontal(
               (pl.col("bp_max") - pl.col("bp_min")) / pl.col("bp_min") * 100,
               (pl.col("rp_max") - pl.col("rp_min")) / pl.col("rp_min") * 100))
           .filter(pl.col("percentage_change") > 100))
    out.write_csv(outpath, include_header=False)
    return out.height


def _run_c(paths, outpath):
    """Decompress + scan in C; marshal the 5 result columns back with bulk copies."""
    lib = ctypes.CDLL(_SO)
    lib.gaia_scan.restype = ctypes.c_long
    lib.gaia_scan.argtypes = [
        ctypes.POINTER(ctypes.c_char_p), ctypes.c_int,
        ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.c_long]
    m = 1 << 18  # upper bound on rows (dataset has ~75k)
    arr = (ctypes.c_char_p * len(paths))(*[p.encode() for p in paths])
    ids = (ctypes.c_int64 * m)()
    bmn = (ctypes.c_double * m)(); bmx = (ctypes.c_double * m)()
    rmn = (ctypes.c_double * m)(); rmx = (ctypes.c_double * m)()
    n = lib.gaia_scan(arr, len(paths), ids, bmn, bmx, rmn, rmx, m)

    def col(buf, typ):
        a = array.array(typ)
        a.frombytes(memoryview(buf).cast("B")[:n * 8])
        return a

    df = pl.DataFrame({"source_id": col(ids, "q"),
                       "bp_min": col(bmn, "d"), "bp_max": col(bmx, "d"),
                       "rp_min": col(rmn, "d"), "rp_max": col(rmx, "d")})
    return _finalize(df, outpath)


def _agg(col):
    """Parse a quoted '[a,b,NaN,...]' flux array into a List(Float64) (fallback)."""
    return (pl.col(col).str.replace_all("NaN", "").str.strip_chars("[]")
            .str.split(",").cast(pl.List(pl.Float64), strict=False))


def _run_polars(paths, outpath):
    from concurrent.futures import ThreadPoolExecutor

    def read_one(path):
        return (pl.read_csv(path, comment_prefix="#",
                            columns=["source_id", "bp_flux", "rp_flux"])
                .select("source_id",
                        bp_min=_agg("bp_flux").list.min(), bp_max=_agg("bp_flux").list.max(),
                        rp_min=_agg("rp_flux").list.min(), rp_max=_agg("rp_flux").list.max()))
    with ThreadPoolExecutor(max_workers=min(len(paths), os.cpu_count() or 8)) as ex:
        frames = list(ex.map(read_one, paths))
    return _finalize(pl.concat(frames), outpath)


def run(indir="/home/irisowner/dev/data/in",
        outpath="/home/irisowner/dev/data/out/result.csv"):
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    # Largest file first so the long-pole decompression starts at t=0.
    paths = sorted(glob.glob(os.path.join(indir, "*.csv.gz")),
                   key=os.path.getsize, reverse=True)
    if _ensure_so():
        return _run_c(paths, outpath)
    return _run_polars(paths, outpath)
