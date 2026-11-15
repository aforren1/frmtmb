"""Review of lane hygiene: codemeta.json against DESCRIPTION (read only)."""
import json
import re
import subprocess
import sys

root = sys.argv[1]
cm = json.load(open(root + "/codemeta.json", encoding="utf-8"))
print("version", cm.get("version"))
print("keys", sorted(cm.keys()))
desc = {}
key = None
for line in open(root + "/DESCRIPTION", encoding="utf-8"):
    line = line.rstrip("\r\n")
    if line[:1] in (" ", "\t") and key:
        desc[key] += " " + line.strip()
    elif ":" in line:
        key, val = line.split(":", 1)
        desc[key] = val.strip()
print("DESCRIPTION Version", desc["Version"])


def deps(field):
    out = {}
    for p in desc.get(field, "").split(","):
        p = p.strip()
        if not p:
            continue
        m = re.match(r"^([A-Za-z0-9.]+)\s*(\((.*)\))?", p)
        out[m.group(1)] = (m.group(3) or "").replace(" ", "")
    return out


want_req = {**deps("Depends"), **deps("Imports")}
want_sug = deps("Suggests")


def cmdeps(k):
    out = {}
    v = cm.get(k, [])
    if isinstance(v, dict):
        v = [x for x in v.values() if isinstance(x, dict)]
    for d in v:
        nm = d.get("name") or d.get("identifier")
        out[nm] = (d.get("version") or "").replace(" ", "")
    return out


have_req = cmdeps("softwareRequirements")
have_sug = cmdeps("softwareSuggestions")
print("requirements: want", len(want_req), "have", len(have_req))
print("  missing", sorted(set(want_req) - set(have_req)))
print("  extra", sorted(set(have_req) - set(want_req)))
print("  version mismatch", {k: (want_req[k], have_req[k]) for k in
                             set(want_req) & set(have_req)
                             if want_req[k] != have_req[k]})
print("suggestions: want", len(want_sug), "have", len(have_sug))
print("  missing", sorted(set(want_sug) - set(have_sug)))
print("  extra", sorted(set(have_sug) - set(want_sug)))
print("  version mismatch", {k: (want_sug[k], have_sug[k]) for k in
                             set(want_sug) & set(have_sug)
                             if want_sug[k] != have_sug[k]})
for k in ("name", "description", "codeRepository", "issueTracker", "url",
          "fileSize", "dateCreated", "dateModified", "datePublished",
          "runtimePlatform", "license", "relatedLink"):
    v = cm.get(k)
    print(k, ":", (json.dumps(v)[:300] if v is not None else None))
print("citation", json.dumps(cm.get("citation"))[:800])
print("DESCRIPTION Title", desc.get("Title"))
print("DESCRIPTION URL", desc.get("URL"))
