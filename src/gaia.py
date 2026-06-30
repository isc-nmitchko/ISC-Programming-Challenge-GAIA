#            .        *        .           .          *       .
#       *        .         ____      _      ___      _      .        *
#   .       *           / ___|    / \    |_ _|    / \         .
#         .            | |  _    / _ \    | |    / _ \    *        .
#     *        .       | |_| |  / ___ \   | |   / ___ \      .         *
#   .              *    \____| /_/   \_\ |___| /_/   \_\  .        *
#        .      *        .          .         *        .        .
#  *   . D R 3   E P O C H   P H O T O M E T R Y   *   .      *    .
# ----------------------------------------------------------------------
"""Gaia DR3 epoch-photometry processor (IRIS Embedded Python + polars).

For every source_id, find the min/max valid BP and RP flux (ignoring NaN/null),
compute (max-min)/min*100 for each band, keep the larger as percentage_change,
and emit the sources whose percentage_change exceeds 100%.

A faster libdeflate + C-kernel variant is available, opt-in, in gaia_fast.py.
"""
import glob, os
from concurrent.futures import ThreadPoolExecutor
import polars as pl


#  *  .   parse a quoted '[a,b,NaN,...]' flux array into a List(Float64)   .  *
#   . NaN -> '' at the string level so the float cast yields null, which .
#  *  .       list.min/list.max already skip -- no per-row list.eval.      .  *
def _agg(col):
    return (pl.col(col).str.replace_all("NaN", "").str.strip_chars("[]")
            .str.split(",").cast(pl.List(pl.Float64), strict=False))


#  .   *    decompress + parse one file, reduce to per-source scalars     *   .
#   *  . polars decompresses .gz natively in Rust; aggregating per file  .   *
#  .       keeps the cross-file concat tiny (5 scalars, not 2 lists).        .
def _read_one(path):
    return (pl.read_csv(path, comment_prefix="#",
                        columns=["source_id", "bp_flux", "rp_flux"])
            .select("source_id",
                    bp_min=_agg("bp_flux").list.min(), bp_max=_agg("bp_flux").list.max(),
                    rp_min=_agg("rp_flux").list.min(), rp_max=_agg("rp_flux").list.max()))


#         *   .        .      ____  __  __  _  _        .   *        .
#    .        *    .        |  _ \|  \/  || \| |   *        .       *
#       fan out across the constellation of files, then reduce to one sky
#    *        .        .    |_| \_|_|\/|_||_|\_|        .        *    .
def run(indir="/home/irisowner/dev/data/in",
        outpath="/home/irisowner/dev/data/out/result.csv"):
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    # Largest file first so the long pole starts at t=0 and finishes within the pool.
    paths = sorted(glob.glob(os.path.join(indir, "*.csv.gz")),
                   key=os.path.getsize, reverse=True)
    with ThreadPoolExecutor(max_workers=min(len(paths), os.cpu_count() or 8)) as ex:
        frames = list(ex.map(_read_one, paths))
    out = (pl.concat(frames)
           .with_columns(percentage_change=pl.max_horizontal(
               (pl.col("bp_max") - pl.col("bp_min")) / pl.col("bp_min") * 100,
               (pl.col("rp_max") - pl.col("rp_min")) / pl.col("rp_min") * 100))
           .filter(pl.col("percentage_change") > 100))
    out.write_csv(outpath, include_header=False)
    return out.height
#       .        *        .        *        .        *        .        *
#    *      .  "we are made of star-stuff"  -- and bytes of it.  .      *
