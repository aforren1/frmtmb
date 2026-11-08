"""Lane ci, last items: stan-next's job limit, continue-on-error on the
record step, and the Stan cache split into restore and save, so that a
job that fails or times out still keeps what it compiled. Each old text
must occur exactly once; line endings are kept as found."""
import pathlib

RESTORE_OLD = """      - uses: actions/cache@v4
        with:
          path: dev/stan-cache
          key: ${{ steps.stankey.outputs.prefix }}${{ hashFiles('dev/stan-cache-key.txt') }}
          restore-keys: ${{ steps.stankey.outputs.prefix }}
"""
RESTORE_NEW = """      # Restore and save are separate steps, so that a job that fails or
      # hits a time limit still keeps the programs it compiled (the
      # combined action saves only from a successful job). A save from
      # such a job goes under the primary key plus "-partial-<run id>",
      # never under the primary key itself: GitHub never overwrites a
      # key, and a partial cache under the primary key would be an exact
      # hit ever after, so the programs it lacks would be compiled on
      # every run and never saved. The restore-keys prefix matches the
      # partial key and takes the newest match, so the next run starts
      # from it and saves the full set under the primary key.
      - uses: actions/cache/restore@v4
        id: stancache
        with:
          path: dev/stan-cache
          key: ${{ steps.stankey.outputs.prefix }}${{ hashFiles('dev/stan-cache-key.txt') }}
          restore-keys: ${{ steps.stankey.outputs.prefix }}
"""
SAVE = """      # last, on any outcome (see the restore step); an exact hit has
      # nothing new to save
      - uses: actions/cache/save@v4
        if: ${{ always() && steps.stankey.outcome == 'success' && steps.stancache.outputs.cache-hit != 'true' && hashFiles('dev/stan-cache/*.rds') != '' }}
        with:
          path: dev/stan-cache
          key: ${{ steps.stankey.outputs.prefix }}${{ hashFiles('dev/stan-cache-key.txt') }}${{ job.status != 'success' && format('-partial-{0}', github.run_id) || '' }}
"""


def edit(path, pairs, append=None):
    p = pathlib.Path(path)
    b = p.read_bytes()
    crlf = b"\r\n" in b
    t = b.decode("utf-8").replace("\r\n", "\n")
    for old, new in pairs:
        assert t.count(old) == 1, (path, old[:60])
        t = t.replace(old, new)
    if append:
        assert t.endswith("\n")
        t = t + append
    if crlf:
        t = t.replace("\n", "\r\n")
    p.write_bytes(t.encode("utf-8"))
    print(path, "crlf" if crlf else "lf")


edit(".github/workflows/brms-likelihood.yaml", [
    ("""      # All five, not two. actions/cache saves at post-job only when the
      # primary key was not an exact hit.""",
     """      # All five, not two. The save step below saves only when the
      # primary key was not an exact hit."""),
    (RESTORE_OLD, RESTORE_NEW),
], append=SAVE)

edit(".github/workflows/ubuntu-next.yaml", [
    ("""      - name: Record the numerical libraries
        env:
""",
     """      # a record only: its failure must not skip the check and the gate
      - name: Record the numerical libraries
        continue-on-error: true
        env:
"""),
    ("""    # Wider than the push job's 90 minutes, because a weekly run finds
    # its cache evicted (GitHub drops a cache unused for 7 days) and
    # compiles every program: all six jobs of 2026-10-19, 10-26 and 11-02
    # missed it. The identity tier alone took 2936, 3826 and 5296 s on
    # 24.04 and 4692, 4688 and 3956 s on 26.04, and at 5296 s the 24.04
    # column was cancelled at 90 minutes
    timeout-minutes: 150
""",
     """    # The push job's 180 minutes, for its reasons, and because a weekly
    # run finds its cache evicted (GitHub drops a cache unused for 7
    # days) and compiles every program: all six jobs of 2026-10-19, 10-26
    # and 11-02 missed it. The identity tier alone took 2936, 3826 and
    # 5296 s on 24.04 and 4692, 4688 and 3956 s on 26.04, and at 5296 s
    # the 24.04 column was cancelled at the old 90 minutes. With a slow
    # mirror's 41 minutes of setup, the tier's 120-minute step limit must
    # fit inside the job's, or a wedged tier cancels the job and skips
    # the tiers after it
    timeout-minutes: 180
"""),
    ("""  # together. Differences from it: the matrix, the version record, the
  # image in the compiled-model cache key, and the job's time limit.
""",
     """  # together. Differences from it: the matrix, the version record, and
  # the image in the compiled-model cache key.
"""),
    (RESTORE_OLD, RESTORE_NEW),
], append=SAVE)
