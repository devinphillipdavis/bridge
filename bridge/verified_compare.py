"""bridge.verified_compare — the budget-matched comparison, adjudicated by osiris_governance.

Requires the optional ``[governance]`` extra. Every number the decision uses is
recomputed from raw Aer counts:

- each side's verification score is re-derived from its verification counts by
  ``ExternalVerifier`` (mean per-qubit |+> survival, the same metric as
  ``dnalang.metrics.survival_plus``) and must match what ``compare`` reported;
- held-out status is checked, not asserted: no verification count set may share a
  ``counts_sha`` with any evaluation the search runs logged in their ledgers;
- the promotion gate's batches (one per seed) come from those packages, and its
  pseudo-replication floor stays at the gate default (5 batches), not the seed count.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set

from osiris_governance.contracts import PermitAuthority
from osiris_governance.evidence_verifier import ExternalVerifier, RawDataPackage
from osiris_governance.ledger import DynamicEvidenceLedger
from osiris_governance.models.ledger import EvidencePlane
from osiris_governance.stat_gate import StatisticalPromotionGate

from dnalang.ledger import Ledger

from .compare import HARD, compare
from .noise import Calibration

METRIC = "mean_qubit_zero_probability"


def counts_digest(counts: Dict[str, int]) -> str:
    """Same digest QuantumFitness logs as ``counts_sha`` for every search evaluation."""
    return hashlib.sha256(json.dumps(counts, sort_keys=True).encode()).hexdigest()[:16]


def search_counts_digests(workdir: Path, seeds: Sequence[int]) -> Set[str]:
    """Every count set the organism and GA searches evaluated, from their hash-chained ledgers."""
    out: Set[str] = set()
    for sd in seeds:
        for name in (f"organism_seed{sd}", f"ga_seed{sd}"):
            led = Ledger(workdir / f"{name}.ledger.jsonl")
            problem = led.verify()
            if problem:
                raise ValueError(f"search ledger {name} is broken: {problem}")
            out.update(e["counts_sha"] for e in led if e.get("kind") == "eval")
    return out


def _package(method: str, row: Dict[str, Any], cal: Calibration, captured: str) -> RawDataPackage:
    side = row[method]
    return RawDataPackage(
        raw_id=f"BRIDGE-{method.upper()}-seed{row['seed']}-{cal.hash()[:8]}",
        source_plane=EvidencePlane.SCIENTIFIC_EXPERIMENTAL,
        provider_receipt_id=None,
        raw_payload={"counts": side["verify_counts"]},
        captured_at_utc=captured,
        metadata={"metric": METRIC, "calibration_hash": cal.hash(), "seed": row["seed"],
                  "genome_key": side["genome_key"]},
    )


def verified_compare(
    workdir: Path,
    seeds: Sequence[int] = range(5),
    budget: int = 300,
    cal: Calibration = HARD,
    tol: float = 0.005,
    verify_shots: int = 1024,
    ledger: Optional[DynamicEvidenceLedger] = None,
    permit_authority: Optional[PermitAuthority] = None,
    alpha: float = 0.05,
) -> Dict[str, Any]:
    """Runs ``compare`` and adjudicates organism (challenger) vs GA (champion) from raw counts."""
    if ledger is not None and permit_authority is None:
        raise ValueError("recording to a ledger needs a permit_authority")
    workdir = Path(workdir)
    res = compare(workdir=workdir, seeds=seeds, budget=budget, cal=cal, tol=tol,
                  verify_shots=verify_shots, keep_counts=True)
    captured = datetime.now(timezone.utc).isoformat()
    verifier = ExternalVerifier(tolerance=1e-9, permit_authority=permit_authority)

    packs: Dict[str, List[RawDataPackage]] = {"organism": [], "ga": []}
    integrity: List[Dict[str, Any]] = []
    for row in res["rows"]:
        for method in packs:
            pkg = _package(method, row, cal, captured)
            sealed = verifier.verify_and_seal(
                raw_pkg=pkg,
                claimed_metrics={METRIC: float(row[method]["verify"])},
                claim_id=f"CLM-{pkg.raw_id}",
                subject=f"{method} verification score, seed {row['seed']}",
            )
            integrity.append({"raw_id": pkg.raw_id, "matched": sealed.all_matched,
                              "recomputed": sealed.recomputed_metrics.get(METRIC)})
            if ledger is not None:
                verifier.record_to_ledger(
                    ledger=ledger, verified_pack=sealed, claim_id=f"CLM-{pkg.raw_id}",
                    subject=f"{method} verification score, seed {row['seed']}",
                    statement=f"{method} verification survival recomputed from raw Aer counts",
                )
            packs[method].append(pkg)

    all_matched = all(i["matched"] for i in integrity)
    search = search_counts_digests(workdir, seeds)
    heldout_pkgs = packs["ga"] + packs["organism"]
    gate = StatisticalPromotionGate(alpha=alpha, mde=0.0, n_resamples=1000)
    decision = gate.evaluate_promotion(
        champion_batches=[verifier.batch_values(p) for p in packs["ga"]],
        challenger_batches=[verifier.batch_values(p) for p in packs["organism"]],
        search_batch_digests=search,
        heldout_batch_digests=[counts_digest(p.raw_payload["counts"]) for p in heldout_pkgs],
    )
    promoted = decision.promoted and all_matched

    out = {
        **res,
        "statistical_promotion_decision": {
            "promoted": promoted,
            "delta_median": decision.delta_median,
            "ci_lower": decision.ci_lower,
            "ci_upper": decision.ci_upper,
            "reason_code": decision.reason_code if all_matched else "EVIDENCE_INTEGRITY_MISMATCH",
            "message": decision.message,
            "search_evaluations_checked": len(search),
        },
        "evidence_integrity": integrity,
        "verified_evidence_sealed": all_matched,
    }
    (workdir / "verified_compare.json").write_text(json.dumps(out, indent=2))
    return out


__all__ = ["verified_compare", "counts_digest", "search_counts_digests"]
