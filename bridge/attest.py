"""
bridge.attest — retrospective attestation of finished drift-check runs
======================================================================

The ledgers of E2–E4 are hash-chained but were not signed while the jobs ran (see "Deviations and
disclosures" in ``results/drift_checks/README.md``). ``build`` writes a manifest that binds, for
each finished run:
- the pre-registration: SHA-256 of the file and the path of its OpenTimestamps proof. The proof is
  not hashed, because ``ots upgrade`` adds attestations to it; it verifies against the file's hash.
- the ledger: SHA-256 of the file, row count, first and last timestamps, and the final chain hash.
  ``build`` refuses a chain that does not verify.
- the telemetry run: run id, record count and head, from ``osiris_governance.telemetry.verify_run``.
- every file in the run's folder, by SHA-256.

An enrolled node key signs the manifest's exact bytes (``osiris_governance.signing``). The signature
is retrospective: it shows that the key holder vouched for these bytes when signing. It does not
show that the ledgers were signed while the jobs ran, that no row is missing, or anything about
what the data mean. ``verify`` checks the signature and recomputes every hash, so a changed,
missing or extra file fails it.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from osiris_governance import signing, telemetry

from bridge.drift_check import sha256_file
from dnalang.ledger import Ledger

DOCUMENT_TYPE = "OSIRIS-RETROSPECTIVE-ATTESTATION-V1"
# run name -> (pre-registration, run folder), both relative to results/drift_checks
RUNS = {"E2": ("PREREGISTRATION.md", "hardware"),
        "E3": ("AMPLITUDE_PREREGISTRATION.md", "amplitude"),
        "E4": ("E4_PREREGISTRATION.md", "pulse_amplitude"),
        "E5": ("E5_PREREGISTRATION.md", "gain_drift")}
STATEMENT = ("Retrospective attestation. The signer vouches that the files listed here are the records "
             "of these runs as committed. The ledgers were not signed while the jobs ran.")
DOES_NOT_SHOW = ["that any ledger row was signed while its job ran",
                 "that no ledger row is missing",
                 "anything about what the data mean"]


def run_entry(root: Path, name: str) -> Dict[str, Any]:
    prereg, folder = RUNS[name]
    path = root / folder / "ledger.jsonl"
    if not path.is_file():                     # Ledger() would create an empty file
        raise FileNotFoundError(f"{name}: no ledger at {path}")
    ledger = Ledger(path)
    problem = ledger.verify()
    if problem:
        raise ValueError(f"{name}: ledger does not verify: {problem}")
    rows = list(ledger)
    if not rows:
        raise ValueError(f"{name}: ledger is empty")
    summary = telemetry.verify_run((root / folder / "telemetry.jsonl").read_bytes().splitlines())
    if not summary.ok:
        raise ValueError(f"{name}: telemetry does not verify: {summary.violations}")
    files = sorted(p for p in (root / folder).rglob("*") if p.is_file())
    return {"preregistration": {"path": prereg, "sha256": sha256_file(root / prereg), "ots": prereg + ".ots"},
            "ledger": {"path": f"{folder}/ledger.jsonl", "sha256": sha256_file(path), "rows": len(rows),
                       "first_ts": rows[0]["ts"], "last_ts": rows[-1]["ts"], "head": rows[-1]["hash"]},
            "telemetry": {"path": f"{folder}/telemetry.jsonl", "run_id": summary.run_id,
                          "records": summary.records, "head": summary.head},
            "files": {p.relative_to(root).as_posix(): sha256_file(p) for p in files}}


def build(out: Path, names: Sequence[str]) -> Dict[str, Any]:
    """Writes the manifest for the named runs to ``out``; paths in it are relative to out's folder."""
    if out.exists():
        raise FileExistsError(f"{out} exists; an attestation is never rewritten")
    manifest = {"document_type": DOCUMENT_TYPE, "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "statement": STATEMENT, "does_not_show": DOES_NOT_SHOW,
                "runs": {name: run_entry(out.parent, name) for name in names}}
    out.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    return manifest


def check(manifest_path: Path) -> List[str]:
    """Recomputes every entry; returns what no longer matches (empty if everything does)."""
    try:
        runs = json.loads(manifest_path.read_text())["runs"]
    except (ValueError, KeyError):
        return ["the manifest is not a readable attestation"]
    problems = []
    for name, recorded in runs.items():
        try:
            now = run_entry(manifest_path.parent, name)
        except (OSError, ValueError) as exc:
            problems.append(str(exc))
            continue
        for part in ("preregistration", "ledger", "telemetry"):
            if now[part] != recorded[part]:
                problems.append(f"{name}: {part} differs")
        old, new = recorded["files"], now["files"]
        problems += [f"{name}: missing {p}" for p in sorted(set(old) - set(new))]
        problems += [f"{name}: extra {p}" for p in sorted(set(new) - set(old))]
        problems += [f"{name}: changed {p}" for p in sorted(set(old) & set(new)) if old[p] != new[p]]
    return problems


def signature_path(manifest_path: Path) -> Path:
    return manifest_path.with_name(manifest_path.name + ".sig")


def sign(manifest_path: Path, private_key_pem: bytes) -> Dict[str, Any]:
    """Signs the manifest's exact bytes and writes ``<manifest>.sig``. Refuses if anything changed."""
    problems = check(manifest_path)
    if problems:
        raise ValueError(f"not signing: {problems}")
    out = signature_path(manifest_path)
    if out.exists():
        raise FileExistsError(f"{out} exists")
    key = signing.load_private_key_pem(private_key_pem)
    data = manifest_path.read_bytes()
    record = {"algorithm": signing.SIGNATURE_ALGORITHM, "manifest": manifest_path.name,
              "manifest_sha256": sha256_file(manifest_path),
              "signer_fingerprint": signing.public_key_fingerprint(key.public_key()),
              "signed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "signature": signing.encode_signature(signing.sign_bytes(key, data))}
    out.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n")
    return record


def verify(manifest_path: Path, public_key_pem: bytes) -> List[str]:
    """Problems with the signature or the files; empty if the signature is valid and nothing changed."""
    problems = []
    record = json.loads(signature_path(manifest_path).read_text())
    key = signing.load_public_key_pem(public_key_pem)
    if record.get("signer_fingerprint") != signing.public_key_fingerprint(key):
        problems.append("the signature names a different key")
    try:
        valid = signing.verify_signature(key, signing.decode_signature(record["signature"]), manifest_path.read_bytes())
    except (KeyError, ValueError):
        valid = False
    if not valid:
        problems.append("the signature does not match the manifest")
    return problems + check(manifest_path)


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(prog="python -m bridge.attest")
    sub = parser.add_subparsers(dest="command", required=True)
    bld = sub.add_parser("build")
    bld.add_argument("--out", type=Path, required=True)
    bld.add_argument("--runs", nargs="+", choices=sorted(RUNS), required=True)
    sgn = sub.add_parser("sign")
    sgn.add_argument("--manifest", type=Path, required=True)
    sgn.add_argument("--key-file", type=Path, required=True, help="the enrolled node key, Ed25519 PEM")
    ver = sub.add_parser("verify")
    ver.add_argument("--manifest", type=Path, required=True)
    ver.add_argument("--public-key", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "build":
        manifest = build(args.out, args.runs)
        print(json.dumps({n: {k: r["ledger"][k] for k in ("rows", "head")} for n, r in manifest["runs"].items()}, indent=1))
    elif args.command == "sign":
        record = sign(args.manifest, args.key_file.read_bytes())
        print(json.dumps({k: record[k] for k in ("manifest", "signer_fingerprint", "signed_at")}, indent=1))
    else:
        problems = verify(args.manifest, args.public_key.read_bytes())
        print(json.dumps({"valid": not problems, "problems": problems}, indent=1))
        raise SystemExit(1 if problems else 0)


if __name__ == "__main__":
    main()
