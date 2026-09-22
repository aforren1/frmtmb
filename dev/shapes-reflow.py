# Re-wrap R lines the shape rewrite pushed past 80 columns.
#
# The house rule is 80 columns in R sources, and the mechanical edits of
# items 2.6d and 2.6f added up to fourteen characters to a call
# (`[, "Estimate"]`) or eight (`_by_dpar`). A break is taken only at a
# comma that is OUTSIDE every string and at the depth of the call being
# broken, and the continuation is indented under the open paren, which
# is what the surrounding code does. A line with no such comma is left
# alone and reported, so that the remainder is a short hand list rather
# than a silent skip.
#
#   python dev/shapes-reflow.py [--apply] <file> ...

import io
import sys


def scan(line):
    """(depth, in_string) after each character, and the open-paren cols."""
    depth = 0
    quote = None
    esc = False
    opens = []          # column of each unclosed open paren
    marks = []          # (index, depth) of every comma outside a string
    for i, ch in enumerate(line):
        if quote is not None:
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == quote:
                quote = None
            continue
        if ch in "\"'":
            quote = ch
        elif ch == "#":
            break        # a comment: nothing after it is code
        elif ch in "([{":
            depth += 1
            opens.append(i)
        elif ch in ")]}":
            depth -= 1
            if opens:
                opens.pop()
        elif ch == ",":
            marks.append((i, depth))
    return marks, opens


def rewrap(line, width=80):
    if len(line) <= width or line.lstrip().startswith("#"):
        return None
    marks, _ = scan(line)
    if not marks:
        return None
    # the shallowest depth a comma sits at is the call being broken
    best = None
    for i, d in marks:
        head = line[:i + 1]
        if len(head) > width:
            continue
        # the indent is the column just after the open paren enclosing
        # this comma
        m2, opens = scan(head)
        if not opens:
            continue
        ind = opens[-1] + 1
        tail = line[i + 1:].lstrip()
        if ind + len(tail) > width:
            continue
        if best is None or i > best[0]:
            best = (i, head.rstrip(), " " * ind + tail)
    if best is None:
        return None
    return [best[1], best[2]]


def main(argv):
    apply = "--apply" in argv
    files = [a for a in argv if not a.startswith("--")]
    left = 0
    done = 0
    for p in files:
        lines = io.open(p, encoding="utf-8").read().split("\n")
        out = []
        changed = False
        for ln in lines:
            new = rewrap(ln)
            if new is None:
                if len(ln) > 80 and not ln.lstrip().startswith("#"):
                    left += 1
                    print("LEFT", p, ln[:70])
                out.append(ln)
            else:
                out.extend(new)
                done += 1
                changed = True
        if changed and apply:
            io.open(p, "w", encoding="utf-8",
                    newline="").write("\n".join(out))
    print("rewrapped", done, "left", left)


if __name__ == "__main__":
    main(sys.argv[1:])
